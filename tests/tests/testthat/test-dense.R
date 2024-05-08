# This tests the dense matrix extraction.
# library(testthat); source("setup.R"); source("test-dense.R")

#################################
# Simple matrices with no chunks. 

set.seed(100000)
mat <- matrix(runif(10000), ncol=200)
big_test_suite(mat, cache.fraction = 0)
