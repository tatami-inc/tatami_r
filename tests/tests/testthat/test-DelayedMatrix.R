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
