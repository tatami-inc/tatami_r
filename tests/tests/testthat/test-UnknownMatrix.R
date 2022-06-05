# Check for correct parsing with SparseArraySeed objects.
# library(testthat); library(raticate.tests); source("test-UnknownMatrix.R")

library(Matrix)
library(DelayedArray)

test_that("Works for UnknownMatrix objects", {
    y <- DelayedArray(matrix(runif(100), 20, 5))
    y <- log1p(y) + rnorm(nrow(y))

    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), 20L)
    expect_identical(raticate.tests::ncol(z), 5L)

    expect_identical(raticate.tests::row(z, 11), y[11,])
    expect_identical(raticate.tests::column(z, 4), y[,4])
})
