# P_value_combination

This is a repository that provides p-value combination functions used in simulations. <br />

R functions include the executation of a python script to calculate Cauchy combination test p-value with equal weights. <br />

To use the function with method="Cauchy", you need to install "reticulate" package in R. This will generate a miniconda environment named "r-reticulate". <br />
Then, 'numpy' and 'mpmath' modules are required. This can be installed in R using following commands. <br />

conda_install("r-reticulate", "numpy") <br />
conda_install("r-reticulate", "mpmath") <br />

The examplary R main script is also provided with virtual p-value matrix. <br />
Users should change the p-value matrix with that generated from their own data. <br />
The matrix should have markers (rows) x methods (columns) format. <br />
When only single marker is provided, covariance of Kost method cannot be calculated. Users should provide the covariance matrix separately. <br />
