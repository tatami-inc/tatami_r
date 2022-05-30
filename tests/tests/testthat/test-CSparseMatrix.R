# Check for correct parsing with CSparseMatrix objects.
# library(testthat); library(raticate.tests); source("test-CSparseMatrix.R")

library(Matrix)

test_that("Works for dgCMatrix objects", {
    x <- rsparsematrix(100, 20, 0.05)

    z <- raticate.tests::parse(x)
    expect_identical(raticate.tests::nrow(z), 100L)
    expect_identical(raticate.tests::ncol(z), 20L)

    expect_identical(raticate.tests::row(z, 5), x[5,])
    expect_identical(raticate.tests::column(z, 9), x[,9])
})

test_that("Works for lgCMatrix objects", {
    x <- rsparsematrix(20, 50, 0.05)

    z <- raticate.tests::parse(x)
    expect_identical(raticate.tests::nrow(z), 20L)
    expect_identical(raticate.tests::ncol(z), 50L)

    expect_identical(raticate.tests::row(z, 12), x[12,])
    expect_identical(raticate.tests::column(z, 29), x[,29])
})
