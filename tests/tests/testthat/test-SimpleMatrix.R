# Check for correct parsing with simple matrix objects.
# library(testthat); library(raticate.tests); source("test-SimpleMatrix.R")

library(Matrix)

test_that("Works for numeric matrix objects", {
    x <- matrix(runif(100), 25, 4)

    z <- raticate.tests::parse(x)
    expect_identical(raticate.tests::nrow(z), 25L)
    expect_identical(raticate.tests::ncol(z), 4L)

    expect_identical(raticate.tests::row(z, 15), x[15,])
    expect_identical(raticate.tests::column(z, 3), x[,3])
})

test_that("Works for integer matrix objects", {
    x <- matrix(rpois(100, 2), 10, 10)

    z <- raticate.tests::parse(x)
    expect_identical(raticate.tests::nrow(z), 10L)
    expect_identical(raticate.tests::ncol(z), 10L)

    expect_identical(raticate.tests::row(z, 2), as.double(x[2,]))
    expect_identical(raticate.tests::column(z, 9), as.double(x[,9]))
})

test_that("Works for logical matrix objects", {
    x <- matrix(rbinom(100, 1, 0.5) > 0, 5, 20)

    z <- raticate.tests::parse(x)
    expect_identical(raticate.tests::nrow(z), 5L)
    expect_identical(raticate.tests::ncol(z), 20L)

    expect_identical(raticate.tests::row(z, 4), as.double(x[4,]))
    expect_identical(raticate.tests::column(z, 19), as.double(x[,19]))
})
