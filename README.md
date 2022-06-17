# P_value_combination

This is a repository that provides p-value combination functions used in simulations 

R functions include the executation of a python script to calculate Cauchy combination test p-value with equal weights. 

To use the function with method="Cauchy", you need to install "reticulate" package in R. This will generate a miniconda environment named "r-reticulate". 
Then, 'numpy' and 'mpmath' modules are required. This can be installed in R using following commands 

conda_install("r-reticulate", "numpy")
conda_install("r-reticulate", "mpmath")

The examplary R main script is also provided with virtual p-value matrix. 
Users should change the p-value matrix with that generated from their own data. 
The matrix should have markers (rows) x methods (columns) format. 
When only single marker is provided, covariance of Kost method cannot be calculated. Users should provide the covariance matrix separately. 
