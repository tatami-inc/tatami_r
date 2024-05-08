# This tests the dense matrix extraction.
# library(testthat); source("setup.R"); source("test-dense-unchunked.R")

set.seed(100000)

{
    NR <- 49
    NC <- 103
    mat <- matrix(runif(NR * NC), ncol=NC)

    test_that("dense unchunked double matrix passes basic checks", {
        expect_type(mat, "double")
        expect_null(DelayedArray::chunkGrid(mat))

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_identical(NR, raticate.tests::nrow(parsed))
        expect_identical(NC, raticate.tests::ncol(parsed))
        expect_false(raticate.tests::is_sparse(parsed))
        expect_false(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, cache.fraction = 0)
    big_test_suite(mat, cache.fraction = 0.01)
    big_test_suite(mat, cache.fraction = 0.1)
}

{
    NR <- 122
    NC <- 53
    mat <- matrix(rpois(NR * NC, lambda=10), ncol=NC)

    test_that("dense unchunked integer matrix passes basic checks", {
        expect_type(mat, "integer")

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_identical(NR, raticate.tests::nrow(parsed))
        expect_identical(NC, raticate.tests::ncol(parsed))
        expect_false(raticate.tests::is_sparse(parsed))
        expect_false(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, cache.fraction = 0)
    big_test_suite(mat, cache.fraction = 0.01)
    big_test_suite(mat, cache.fraction = 0.1)
}

{
    NR <- 212
    NC <- 23
    mat <- matrix(rbinom(NR * NC, 1, 0.5) == 1, ncol=NC)

    test_that("dense unchunked logical matrix passes basic checks", {
        expect_type(mat, "logical")

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_identical(NR, raticate.tests::nrow(parsed))
        expect_identical(NC, raticate.tests::ncol(parsed))
        expect_false(raticate.tests::is_sparse(parsed))
        expect_false(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, cache.fraction = 0)
    big_test_suite(mat, cache.fraction = 0.01)
    big_test_suite(mat, cache.fraction = 0.1)
}
