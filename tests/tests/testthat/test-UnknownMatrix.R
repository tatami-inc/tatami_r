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
    expect_identical(raticate.tests::row_subset(z, 10, 2, 4), y[10,2:4])
    expect_identical(raticate.tests::column_subset(z, 3, 11, 15), y[11:15,3])

    expect_identical(raticate.tests::rows(z), lapply(seq_len(nrow(y)), function(i) y[i,]))
    expect_identical(raticate.tests::rows_subset(z, 3, 5), lapply(seq_len(nrow(y)), function(i) y[i,3:5]))
    expect_identical(raticate.tests::columns(z), lapply(seq_len(ncol(y)), function(i) y[,i]))
    expect_identical(raticate.tests::columns_subset(z, 5, 16), lapply(seq_len(ncol(y)), function(i) y[5:16,i]))
})
