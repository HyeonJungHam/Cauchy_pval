# P_value_combination

This is a repository that provides p-value combination functions used in simulations. <br />

**Please put all files in the same directory and execute "p_value_combination_function_main.R" script. <br />**

R functions include the execution of a python script to calculate Cauchy combination test p-value with equal weights. <br />

To use the function with method="Cauchy", you need to install "reticulate" package in R. This will generate a miniconda environment named "r-reticulate". <br />

Then, 'numpy' and 'mpmath' modules are required. This can be installed in R using following commands. <br />

<pre>
library(reticulate) 
conda_install("r-reticulate", "numpy") 
conda_install("r-reticulate", "mpmath")
</pre>

"p_value_combination_function_main.R" includes a virtual p-value matrix. <br />

Users should change the p-value matrix with that generated from their own data. <br />

The matrix should have markers (rows) x methods (columns) format. <br />

When only single marker is provided, covariance of Kost method cannot be calculated. Users should provide the covariance matrix separately. <br />

# **Usage**
<pre>
combine_pvalue <- function(mat, 
                           method = c("Fisher", "MinP", "Simes", "Stouffer", "Kost", "Cauchy"), 
                           cov = NULL, 
                           Cauchy_input_path = "./Cauchy_input/", 
                           Cauchy_output_path = "./Cauchy_output/")
</pre>
