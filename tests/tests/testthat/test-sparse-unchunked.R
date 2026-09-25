# library(testthat); source("setup.R"); source("test-sparse-unchunked.R")

set.seed(100000)

{
    NR <- 34
    NC <- 87
    mat <- as(Matrix::rsparsematrix(NR, NC, 0.15), "SVT_SparseMatrix")
    name <- "sparse unchunked double matrix"

    test_that(paste(name, "passes basic checks"), {
        expect_identical(DelayedArray::type(mat), "double")
        expect_null(DelayedArray::chunkGrid(mat))

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_true(raticate.tests::sparse(parsed))
        expect_false(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, name)
}

{
    NR <- 67
    NC <- 24 
    mat <- matrix(0L, NR, NC)
    nnz <- length(mat) * 0.2
    mat[sample(length(mat), nnz)] <- rpois(nnz, lambda=10)
    mat <- as(mat, "SVT_SparseMatrix")
    name <- "sparse unchunked integer matrix"

    test_that(paste(name, "passes basic checks"), {
        expect_identical(DelayedArray::type(mat), "integer")
        expect_null(DelayedArray::chunkGrid(mat))

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_true(raticate.tests::sparse(parsed))
        expect_false(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, name)
}

# This implicitly checks the all-1 lacunar leaf nodes, where nzvals is set to NULL.
{
    NR <- 151
    NC <- 7
    mat <- as(matrix(rbinom(NR * NC, 1, 0.2) == 1, ncol=NC), "SVT_SparseMatrix")
    name <- "sparse unchunked logical matrix"

    test_that(paste(name, "passes basic checks"), {
        expect_identical(DelayedArray::type(mat), "logical")
        expect_null(DelayedArray::chunkGrid(mat))

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_true(raticate.tests::sparse(parsed))
        expect_false(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, name)
}

{
    NR <- 26
    NC <- 138
    mat <- matrix(0, NR, NC)
    mat[,1:NC %% 2 == 1] <- matrix(rpois(NR * NC / 2, lambda = 2), nrow = NR)
    mat <- as(mat, "SVT_SparseMatrix")
    name <- "sparse unchunked double matrix with empty columns"

    test_that(paste(name, "passes basic checks"), {
        expect_identical(DelayedArray::type(mat), "double")
        expect_null(DelayedArray::chunkGrid(mat))

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_true(raticate.tests::sparse(parsed))
        expect_false(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, name)
}

{
    NR <- 56
    NC <- 87 
    mat <- matrix(0L, NR, NC)
    mat <- as(mat, "SVT_SparseMatrix")
    name <- "sparse unchunked integer matrix with no values"

    test_that(paste(name, "passes basic checks"), {
        expect_identical(DelayedArray::type(mat), "integer")
        expect_null(DelayedArray::chunkGrid(mat))

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_true(raticate.tests::sparse(parsed))
        expect_false(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, name)
}
