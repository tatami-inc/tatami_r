# This tests the dense matrix extraction.
# library(testthat); source("setup.R"); source("test-sparse-unchunked.R")

set.seed(100000)

{
    NR <- 67
    NC <- 173
    mat <- as(Matrix::rsparsematrix(NR, NC, 0.15), "SVT_SparseMatrix")
    expect_identical(DelayedArray::type(mat), "double")
    expect_null(DelayedArray::chunkGrid(mat))

    parsed <- raticate.tests::parse(mat, 0, FALSE)
    expect_true(raticate.tests::is_sparse(parsed))
    expect_false(raticate.tests::prefer_rows(parsed))

    big_test_suite(mat, cache.fraction = 0)
    big_test_suite(mat, cache.fraction = 0.01)
    big_test_suite(mat, cache.fraction = 0.1)
}

{
    NR <- 133
    NC <- 49
    mat <- matrix(0L, NR, NC)
    nnz <- length(mat) * 0.2
    mat[sample(length(mat), nnz)] <- rpois(nnz, lambda=10)
    mat <- as(mat, "SVT_SparseMatrix")
    expect_identical(DelayedArray::type(mat), "integer")

    parsed <- raticate.tests::parse(mat, 0, FALSE)
    expect_true(raticate.tests::is_sparse(parsed))
    expect_false(raticate.tests::prefer_rows(parsed))

    big_test_suite(mat, cache.fraction = 0)
    big_test_suite(mat, cache.fraction = 0.01)
    big_test_suite(mat, cache.fraction = 0.1)
}

{
    NR <- 302
    NC <- 13
    mat <- as(matrix(rbinom(NR * NC, 1, 0.2) == 1, ncol=NC), "SVT_SparseMatrix")
    expect_identical(DelayedArray::type(mat), "logical")

    big_test_suite(mat, cache.fraction = 0)
    big_test_suite(mat, cache.fraction = 0.01)
    big_test_suite(mat, cache.fraction = 0.1)
}
