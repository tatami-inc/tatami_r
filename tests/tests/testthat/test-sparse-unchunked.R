# This tests the dense matrix extraction.
# library(testthat); source("setup.R"); source("test-sparse-unchunked.R")

set.seed(100000)

{
    NR <- 34
    NC <- 87
    mat <- as(Matrix::rsparsematrix(NR, NC, 0.15), "SVT_SparseMatrix")

    test_that("sparse unchunked double matrix passes basic checks", {
        expect_identical(DelayedArray::type(mat), "double")
        expect_null(DelayedArray::chunkGrid(mat))

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_true(raticate.tests::sparse(parsed))
        expect_false(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat)
}

{
    NR <- 67
    NC <- 24 
    mat <- matrix(0L, NR, NC)
    nnz <- length(mat) * 0.2
    mat[sample(length(mat), nnz)] <- rpois(nnz, lambda=10)
    mat <- as(mat, "SVT_SparseMatrix")

    test_that("sparse unchunked integer matrix passes basic checks", {
        expect_identical(DelayedArray::type(mat), "integer")
        expect_null(DelayedArray::chunkGrid(mat))

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_true(raticate.tests::sparse(parsed))
        expect_false(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat)
}

{
    NR <- 151
    NC <- 7
    mat <- as(matrix(rbinom(NR * NC, 1, 0.2) == 1, ncol=NC), "SVT_SparseMatrix")

    test_that("sparse unchunked logical matrix passes basic checks", {
        expect_identical(DelayedArray::type(mat), "logical")
        expect_null(DelayedArray::chunkGrid(mat))

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_true(raticate.tests::sparse(parsed))
        expect_false(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat)
}
