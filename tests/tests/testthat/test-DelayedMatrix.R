# Check for correct parsing with SparseArraySeed objects.
# library(testthat); library(raticate.tests); source("test-DelayedMatrix.R")

library(Matrix)
library(DelayedArray)

test_that("Works for simple DelayedMatrix objects", {
    y <- DelayedArray(rsparsematrix(100, 20, 0.05))

    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), 100L)
    expect_identical(raticate.tests::ncol(z), 20L)

    expect_identical(raticate.tests::row(z, 5), y[5,])
    expect_identical(raticate.tests::column(z, 9), y[,9])
})

test_that("Works for subsetted DelayedMatrix objects", {
    x <- DelayedArray(matrix(runif(300), 15, 20))

    # Only rows.
    y <- x[10:5,]
    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), 6L)
    expect_identical(raticate.tests::ncol(z), 20L)
    expect_identical(raticate.tests::row(z, 4), as.double(y[4,]))
    expect_identical(raticate.tests::column(z, 2), as.double(y[,2]))

    # Only columns.
    y <- x[,c(1,3,5,7,9)]
    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), 15L)
    expect_identical(raticate.tests::ncol(z), 5L)
    expect_identical(raticate.tests::row(z, 8), as.double(y[8,]))
    expect_identical(raticate.tests::column(z, 3), as.double(y[,3]))

    # Both  
    y <- x[3:5, 2:10]
    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), 3L)
    expect_identical(raticate.tests::ncol(z), 9L)
    expect_identical(raticate.tests::row(z, 1), as.double(y[1,]))
    expect_identical(raticate.tests::column(z, 6), as.double(y[,6]))
})

test_that("Works for dimnames-altered objects", {
    y <- DelayedArray(rsparsematrix(100, 20, 0.05))
    dimnames(y) <- list(1:100, LETTERS[1:20])

    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), 100L)
    expect_identical(raticate.tests::ncol(z), 20L)

    expect_identical(raticate.tests::row(z, 5), unname(y[5,]))
    expect_identical(raticate.tests::column(z, 9), unname(y[,9]))
})

test_that("Works for transposed objects", {
    y <- DelayedArray(rsparsematrix(100, 20, 0.05))
    y <- t(y)

    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), 20L)
    expect_identical(raticate.tests::ncol(z), 100L)

    expect_identical(raticate.tests::row(z, 1), unname(y[1,]))
    expect_identical(raticate.tests::column(z, 19), unname(y[,19]))
})

test_that("Works for combined objects", {
    # By row.
    x1 <- DelayedArray(rsparsematrix(50, 20, 0.05))
    x2 <- DelayedArray(matrix(runif(200), 10, 20))
    y <- BiocGenerics::rbind(x1, x2)

    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), 60L)
    expect_identical(raticate.tests::ncol(z), 20L)

    expect_identical(raticate.tests::row(z, 3), unname(y[3,]))
    expect_identical(raticate.tests::row(z, 55), unname(y[55,]))
    expect_identical(raticate.tests::column(z, 16), unname(y[,16]))

    # By row.
    x1 <- DelayedArray(matrix(runif(150), 10, 15))
    x2 <- DelayedArray(rsparsematrix(10, 20, 0.05))
    y <- BiocGenerics::cbind(x1, x2)

    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), 10L)
    expect_identical(raticate.tests::ncol(z), 35L)

    expect_identical(raticate.tests::column(z, 5), unname(y[,5]))
    expect_identical(raticate.tests::column(z, 30), unname(y[,30]))
    expect_identical(raticate.tests::row(z, 8), unname(y[8,]))
})
