# This tests the dense matrix extraction.
# library(testthat); source("setup.R"); source("test-matrix.R")

set.seed(100000)
NR <- 49
NC <- 103
mat <- matrix(runif(NR * NC), ncol=NC)
expect_type(mat, "double")
big_test_suite(mat, cache.fraction = 0)
big_test_suite(mat, cache.fraction = 0.01)
big_test_suite(mat, cache.fraction = 0.1)

NR <- 122
NC <- 53
mat <- matrix(rpois(NR * NC, lambda=10), ncol=NC)
expect_type(mat, "integer")
big_test_suite(mat, cache.fraction = 0)
big_test_suite(mat, cache.fraction = 0.01)
big_test_suite(mat, cache.fraction = 0.1)

NR <- 212
NC <- 23
mat <- matrix(rbinom(NR * NC, 1, 0.5) == 1, ncol=NC)
expect_type(mat, "logical")
big_test_suite(mat, cache.fraction = 0)
big_test_suite(mat, cache.fraction = 0.01)
big_test_suite(mat, cache.fraction = 0.1)
