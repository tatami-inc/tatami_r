# This tests the dense matrix extraction.
# library(testthat); source("setup.R"); source("test-dense-unchunked.R")

set.seed(100000)

{
    NR <- 23
    NC <- 52
    mat <- matrix(runif(NR * NC), ncol=NC)

    test_that("dense unchunked double matrix passes basic checks", {
        expect_type(mat, "double")
        expect_null(DelayedArray::chunkGrid(mat))

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_equal(NR, raticate.tests::num_rows(parsed))
        expect_equal(NC, raticate.tests::num_columns(parsed))
        expect_false(raticate.tests::sparse(parsed))
        expect_false(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, cache.fraction = 0)
    big_test_suite(mat, cache.fraction = 0.01)
    big_test_suite(mat, cache.fraction = 0.1)
}

{
    NR <- 61
    NC <- 27
    mat <- matrix(rpois(NR * NC, lambda=10), ncol=NC)

    test_that("dense unchunked integer matrix passes basic checks", {
        expect_type(mat, "integer")
        expect_null(DelayedArray::chunkGrid(mat))

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_equal(NR, raticate.tests::num_rows(parsed))
        expect_equal(NC, raticate.tests::num_columns(parsed))
        expect_false(raticate.tests::sparse(parsed))
        expect_false(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, cache.fraction = 0)
    big_test_suite(mat, cache.fraction = 0.01)
    big_test_suite(mat, cache.fraction = 0.1)
}

{
    NR <- 106
    NC <- 12
    mat <- matrix(rbinom(NR * NC, 1, 0.5) == 1, ncol=NC)

    test_that("dense unchunked logical matrix passes basic checks", {
        expect_type(mat, "logical")
        expect_null(DelayedArray::chunkGrid(mat))

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_equal(NR, raticate.tests::num_rows(parsed))
        expect_equal(NC, raticate.tests::num_columns(parsed))
        expect_false(raticate.tests::sparse(parsed))
        expect_false(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, cache.fraction = 0)
    big_test_suite(mat, cache.fraction = 0.01)
    big_test_suite(mat, cache.fraction = 0.1)
}
