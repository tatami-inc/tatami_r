# Check for correct parsing with SparseArraySeed objects.
# library(testthat); library(raticate.tests); source("test-SparseArraySeed.R")

library(Matrix)
library(DelayedArray)

test_that("Works for numeric SparseArraySeed objects", {
    x <- rsparsematrix(100, 20, 0.05)
    y <- as(x, "SparseArraySeed")

    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), 100L)
    expect_identical(raticate.tests::ncol(z), 20L)

    expect_identical(raticate.tests::row(z, 5), x[5,])
    expect_identical(raticate.tests::column(z, 9), x[,9])
})

test_that("Works for integer SparseArraySeed objects", {
    x <- rsparsematrix(40, 50, 0.05)
    y <- as(x, "SparseArraySeed")
    y@nzdata <- as.integer(ceiling(y@nzdata * 10))

    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), 40L)
    expect_identical(raticate.tests::ncol(z), 50L)

    ref <- DelayedArray(y)
    expect_identical(raticate.tests::row(z, 7), as.double(ref[7,]))
    expect_identical(raticate.tests::column(z, 2), as.double(ref[,2]))
})

test_that("Works for logical SparseArraySeed objects", {
    x <- rsparsematrix(25, 25, 0.05) != 0
    y <- as(x, "SparseArraySeed")

    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), 25L)
    expect_identical(raticate.tests::ncol(z), 25L)

    expect_identical(raticate.tests::row(z, 10), as.double(x[10,]))
    expect_identical(raticate.tests::column(z, 10), as.double(x[,10]))
})
