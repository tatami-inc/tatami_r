# Check for correct handling of UnknownMatrix objects.
# library(testthat); library(raticate.tests); source("setup.R"); source("test-UnknownMatrix.R")

# TODO: add test for default block size.

#test_that("Behaves correctly with R-side errors", {
#    setClass("MyFailMatrix", contains="dgeMatrix")
#    setMethod("extract_array", "MyFailMatrix", function(x, index) stop("HEY!"))
#    setMethod("type", "MyFailMatrix", function(x) "double") # avoid errors in construction
#    y <- new("MyFailMatrix", Matrix(runif(1000), 50, 20))
#
#    z <- raticate.tests::parse(y)
#    expect_identical(raticate.tests::nrow(z), 50L)
#    expect_error(raticate.tests::row(z, 1), "HEY")
#    expect_error(raticate.tests::rows(z), "HEY")
#    expect_error(raticate.tests::rowsums(z)) # error message changes depending on TEST_CUSTOM_PARALLEL, so don't worry about it.
#
#    # Also behaves correctly if we emit incorrect dimensions.
#    setMethod("extract_array", "MyFailMatrix", function(x, index) {
#        stuff <- callNextMethod()
#        matrix(0, 123, 123)
#    })
#
#    expect_error(raticate.tests::row(z, 1), "incorrect dimensions")
#    expect_error(raticate.tests::rows(z), "incorrect dimensions")
#    expect_error(raticate.tests::rowsums(z)) # error message changes depending on TEST_CUSTOM_PARALLEL, so don't worry about it.
#})
