# Check for correct parsing with SVT_SparseMatrix objects.
# library(testthat); library(raticate.tests); source("setup.R"); source("test-SVT_SparseMatrix.R")

library(Matrix)
library(DelayedArray)

test_that("Works for numeric SVT_SparseMatrix objects", {
    x <- rsparsematrix(100, 20, 0.05)
    y <- as(x, "SVT_SparseMatrix")

    # Check that the class is preserved.
    out <- extract_sparse_array(y, list(NULL, NULL))
    expect_s4_class(out, "SVT_SparseMatrix")

    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), 100L)
    expect_identical(raticate.tests::ncol(z), 20L)

    expect_identical(raticate.tests::row(z, 5), x[5,])
    expect_identical(raticate.tests::column(z, 9), x[,9])
})

test_that("Works for integer SVT_SparseMatrix objects", {
    x <- as.matrix(ceiling(rsparsematrix(40, 50, 0.05) * 10))
    storage.mode(x) <- "integer"
    y <- as(x, "SVT_SparseMatrix")
    expect_identical(type(y), "integer")

    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), 40L)
    expect_identical(raticate.tests::ncol(z), 50L)

    ref <- DelayedArray(y)
    expect_identical(raticate.tests::row(z, 7), as.double(ref[7,]))
    expect_identical(raticate.tests::column(z, 2), as.double(ref[,2]))
})

test_that("Works for logical SVT_SparseMatrix objects", {
    x <- rsparsematrix(25, 25, 0.05) != 0
    y <- as(x, "SVT_SparseMatrix")
    expect_identical(type(y), "logical")

    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), 25L)
    expect_identical(raticate.tests::ncol(z), 25L)

    expect_identical(raticate.tests::row(z, 10), as.double(x[10,]))
    expect_identical(raticate.tests::column(z, 10), as.double(x[,10]))
})
