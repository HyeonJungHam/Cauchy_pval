path = "Current directory where 'Cauchy_pval.py' exists"
setwd(path)

source("p_value_combination_function.R")

mat = matrix(runif(25), nrow = 5, ncol = 5)

c("Fisher", "MinP", "Simes", "Stouffer", "Kost", "Cauchy")

combine_pvalue(mat, method = "Fisher")
combine_pvalue(mat, method = "MinP")
combine_pvalue(mat, method = "Simes")
combine_pvalue(mat, method = "Stouffer")
combine_pvalue(mat, method = "Kost")
combine_pvalue(mat, method = "Cauchy")
