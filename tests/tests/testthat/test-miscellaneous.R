# Check for correct handling of UnknownMatrix objects.
# library(testthat); source("setup.R"); source("test-miscellaneous.R")

library(Matrix)
test_that("works correctly with the default cache size", {
    y <- Matrix(runif(1000), 50, 20)
    z <- raticate.tests::parse(y, -1, FALSE)
    expect_identical(raticate.tests::num_rows(z), 50L)
    expect_identical(raticate.tests::num_columns(z), 20L)
})

# parallelization handles edge cases correctly
{
    library(Matrix)
    y <- Matrix(runif(1), 1, 1) # only one job in either dimension.
    parallel_test_suite(y, 0.1)

    y <- Matrix(runif(0), 0, 0) # no jobs in either dimension.
    parallel_test_suite(y, 0.1)

    y <- Matrix(runif(4), 2, 2) # fewer jobs than threads (for workers = 3).
    parallel_test_suite(y, 0.1)
}

test_that("Behaves correctly with R-side errors", {
    setClass("MyFailMatrix", contains="dgeMatrix")
    setMethod("extract_array", "MyFailMatrix", function(x, index) stop("HEY!"))
    setMethod("type", "MyFailMatrix", function(x) "double") # avoid errors in construction
    y <- new("MyFailMatrix", Matrix(runif(1000), 50, 20))

    z <- raticate.tests::parse(y, 0, FALSE)
    expect_identical(raticate.tests::num_rows(z), 50L)
    expect_error(raticate.tests::myopic_dense_full(z, TRUE, 10), "HEY")
    expect_error(raticate.tests::myopic_sparse_full(z, TRUE, 10, TRUE, TRUE), "HEY")
    expect_error(raticate.tests::oracular_dense_full(z, TRUE, 10), "HEY")
    expect_error(raticate.tests::oracular_sparse_full(z, TRUE, 10, TRUE, TRUE), "HEY")
    expect_error(raticate.tests::myopic_dense_sums(z, TRUE, 1)) # error message changes depending on TEST_CUSTOM_PARALLEL, so don't worry about it.
    expect_error(raticate.tests::myopic_sparse_sums(z, TRUE, 3))

    # Also behaves correctly if we emit incorrect dimensions.
    setMethod("extract_array", "MyFailMatrix", function(x, index) {
        stuff <- callNextMethod()
        matrix(0, 123, 123)
    })

    expect_error(raticate.tests::myopic_dense_full(z, TRUE, 10), "incorrect dimensions")
    expect_error(raticate.tests::myopic_sparse_full(z, TRUE, 10, TRUE, TRUE), "incorrect dimensions")
    expect_error(raticate.tests::oracular_dense_sums(z, TRUE, 1)) # error message changes depending on TEST_CUSTOM_PARALLEL, so don't worry about it.
    expect_error(raticate.tests::oracular_sparse_sums(z, TRUE, 3))
})

test_that("executor setting works as expected", {
    expect_true(raticate.tests::test_set_executor())
})
