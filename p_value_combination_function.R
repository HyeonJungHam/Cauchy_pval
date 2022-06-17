# library(reticulate)
# conda_install("r-reticulate", "numpy")
# conda_install("r-reticulate", "mpmath")
# Let columns be methods to combine and rows be markers

combine_pvalue <- function(mat, 
                           method = c("Fisher", "MinP", "Simes", "Stouffer", "Kost", "Cauchy"), 
                           cov = NULL, 
                           Cauchy_input_path = "./Cauchy_input/", 
                           Cauchy_output_path = "./Cauchy_output/")
{
  na.method <- pmatch(method, c("Fisher", "MinP", "Simes", "Stouffer", "Kost", "Cauchy"))
  if(is.na(na.method)) stop("invalid 'method' argument")
  method <- match.arg(method)
  
  if(!dir.exists(Cauchy_input_path)){
    dir.create(Cauchy_input_path)
  }
  
  if(!dir.exists(Cauchy_output_path)){
    dir.create(Cauchy_output_path)
  }
  
  if(method == "Fisher"){
    ret = compute_fisher(mat)
  }
  else if(method == "MinP"){
    ret = compute_minP(mat)
  }
  else if(method == "Kost"){
    if(nrow(mat)>=1){
      ret = compute_kost(mat, mat)
    }
    else{
      if(missing(cov)){
        stop("Data matrix for calculating covariance needs to be provided to compute Kost method")
      }
      ret = compute_kost(mat, cov)
    }
  }
  else if(method == "Simes"){
    ret = compute_simes(mat)
  }
  else if(method == "Stouffer"){
    ret = compute_stouffer(mat)
  }
  else if(method == "Cauchy"){
    # Cauchy_input_path = 
    compute_cauchy_input(mat, input_path = Cauchy_input_path)
    command = paste0("Cauchy_pval.py -i ", Cauchy_input_path, " -o ", Cauchy_output_path)
    system(paste0("python Cauchy_pval.py -i ", Cauchy_input_path, " -o ", Cauchy_output_path))
    Sys.sleep(5) 
    ret = compute_cauchy(output_path = Cauchy_output_path)
  }
  ret
}

fisher_pvalue = function(marker_p_values){
  marker_p_values = data.frame(marker_p_values[!is.na(marker_p_values)])
  new_T = -2*sum(log(marker_p_values), na.rm = T)
  new_p = pchisq(new_T, df=2*nrow(marker_p_values), lower.tail=FALSE)
  if(new_p==0){
    #print("new_p is 0")
    new_p=1e-323
  }
  new_p
}

compute_fisher = function(mat){
  final_p_list <- c()
  for(marker in 1:nrow(mat)){
    marker_p_values <- mat[marker, ]
    
    # fisher's method 
    new_p_list <- c()
    new_p = fisher_pvalue(marker_p_values)
    final_p_list <- c(final_p_list, new_p)
  }
  return(final_p_list)
}

minP_pvalue = function(marker_p_values){
  T_min = min(marker_p_values, na.rm = T)
  
  p_min = pbeta(T_min, 1, length(marker_p_values))
  p_min
}

compute_minP <- function(mat){
  final_p_list <- c()
  
  marker = 1
  for(marker in 1:nrow(mat)){
    marker_p_values <- mat[marker,]
    
    # min p-value method 
    p_min = minP_pvalue(marker_p_values)
    
    final_p_list <- c(final_p_list, p_min)
  }
  
  return(final_p_list)
}

# Kost method adopted from "EmpiricalBrownsMethod" library 
kostPolyFit <- function(cor) {
  a1 <- 3.263
  a2 <- 0.710
  a3 <- 0.027 #Kost cubic coeficients
  (a1*cor + a2*cor^2 + a3*cor^3)
}

combinePValues <- function(covar_matrix, p_values, extra_info = FALSE){
  N = ncol(covar_matrix) # number of samples
  df_fisher = 2.0*N
  Expected  = 2.0*N
  cov_sum <- (2*sum(covar_matrix[lower.tri(covar_matrix, diag=FALSE)], na.rm = T))
  Var = 4.0*N+cov_sum
  c = Var/(2.0*Expected)
  df_brown = (2.0*Expected^2)/Var
  if (df_brown > df_fisher) {
    df_brown = df_fisher
    c = 1.0
  }
  x = 2.0*sum( -log(p_values), na.rm = T)
  p_brown = pchisq(df=df_brown, q=x/c, lower.tail=FALSE)
  p_fisher = pchisq(df=df_fisher, q=x, lower.tail=FALSE)
  if (extra_info){
    return(list(P_test=p_brown, P_Fisher=p_fisher, Scale_Factor_C=c, DF=df_brown))
  }else{
    return(p_brown)
  }
}

calculateKostCovariance <- function(data_matrix){
  m = nrow(data_matrix)
  covar_matrix = mat.or.vec(m, m)
  for (i in 1:m) {
    for (j in i:m) {
      res0 <- cor.test(data_matrix[i,], data_matrix[j,])
      cor <- res0$estimate
      p_val <- res0$p.value
      covar = kostPolyFit(cor)
      covar_matrix[i, j] = covar
      covar_matrix[j, i] = covar
    }
  }
  return(covar_matrix)
}

compute_kost <- function(mat, cov){
  data_matrix = as.matrix(mat)
  covar_matrix <- calculateKostCovariance(t(cov))
  p_kost = combinePValues(covar_matrix, data_matrix[1, ])
  
  final_T_list <- c()
  final_p_list <- c()
  threshold_list <- c()
  for(marker in 1:nrow(mat)){
    p_values = data_matrix[marker, ]
    N = ncol(covar_matrix) # number of samples
    df_fisher = 2.0*N
    Expected  = 2.0*N
    cov_sum <- (2*sum(covar_matrix[lower.tri(covar_matrix, diag=FALSE)], na.rm = T))
    Var = 4.0*N+cov_sum
    c = Var/(2.0*Expected)
    df_brown = (2.0*Expected^2)/Var
    if (df_brown > df_fisher) {
      df_brown = df_fisher
      c = 1.0
    }
    x = 2.0*sum( -log(p_values), na.rm = T)
    p_brown = pchisq(df=df_brown, q=x/c, lower.tail=FALSE)
    p_fisher = pchisq(df=df_fisher, q=x, lower.tail=FALSE)
    final_T_list <- c(final_T_list, x/c)
    final_p_list <- c(final_p_list, p_brown)
    
    threshold_p_values = t(data.frame(rep(0.05, length(p_values))))
    threshold_x = 2.0*sum( -log(threshold_p_values), na.rm = T)
    threshold_p_brown =  pchisq(df=df_brown, q=threshold_x/c, lower.tail=FALSE)
    threshold_list <- c(threshold_list, threshold_p_brown)
  }
  
  return(final_p_list)
}

simes_pvalue = function(marker_p_values){
  marker_p_values = data.frame(marker_p_values[!is.na(marker_p_values)])
  r=rank(marker_p_values)
  pval=min(nrow(marker_p_values)*marker_p_values/r)
}

compute_simes = function(mat){
  final_p_list <- c()
  for(marker in 1:nrow(mat)){
    marker_p_values <- mat[marker, ]
    
    # fisher's method 
    new_p_list <- c()
    
    new_p = simes_pvalue(marker_p_values)
    
    final_p_list <- c(final_p_list, new_p)
    
  }
  return(final_p_list)
}


stouffer_pvalue = function(marker_p_values){
  marker_p_values = (marker_p_values[!is.na(marker_p_values)])
  w = rep(1, length(marker_p_values))
  zp <- qnorm(marker_p_values, lower.tail = F, log.p = F) %*% 
    w / sqrt(sum(w^2))
  p.val = pnorm(zp, lower.tail = F, log.p = F)
  return(p.val)
}

compute_stouffer = function(mat){
  final_p_list <- c()
  for(marker in 1:nrow(mat)){
    marker_p_values <- mat[marker, ]
    
    # fisher's method 
    new_p_list <- c()
    
    new_p = stouffer_pvalue(marker_p_values)
    
    final_p_list <- c(final_p_list, new_p)
    
  }
  return(final_p_list)
}


compute_cauchy_input <- function(mat, input_path = "./cauchy_input/"){
  
  cauchy_path = input_path
  
  if(!dir.exists(cauchy_path)){
    dir.create(cauchy_path)
  }
  data_matrix = as.matrix(mat)
  file_name = paste0(cauchy_path, "data_matrix_Cauchy.csv")
  write.csv(data_matrix, file_name, row.names = F, quote = F)
  
}

compute_cauchy <- function(output_path = "./cauchy_output/"){
  cauchy_path = output_path
  file_name = paste0(cauchy_path, "data_matrix_Cauchy.csv")
  Cauchy_python = read.csv(file_name, header = F, stringsAsFactors = F)
  
  final_T_list = Cauchy_python$V1
  final_p_list = Cauchy_python$V2
  
  return(final_p_list)
}

