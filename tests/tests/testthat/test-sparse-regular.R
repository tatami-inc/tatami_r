# This tests the dense matrix extraction with regular grids.
# library(testthat); source("setup.R"); source("test-sparse-regular.R")

setClass("RegularChunkedSparseMatrix", contains="SVT_SparseMatrix", slots=c(chunks="integer"))
setMethod("chunkdim", "RegularChunkedSparseMatrix", function(x) x@chunks)
RegularChunkedSparseMatrix <- function(mat, chunks) {
    spmat <- as(mat, "SVT_SparseMatrix")
    new("RegularChunkedSparseMatrix", spmat, chunks=as.integer(chunks))
}

set.seed(150000)

{
    NR <- 57
    NC <- 208
    mat <- RegularChunkedSparseMatrix(Matrix::rsparsematrix(NR, NC, 0.24), chunks=c(11, 7)) 

    test_that("sparse regularly-chunked double matrix passes basic checks", {
        expect_true(is_sparse(mat))
        expect_s4_class(chunkGrid(mat), "RegularArrayGrid")

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_true(raticate.tests::is_sparse(parsed))
        expect_false(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, cache.fraction = 0.01)
    big_test_suite(mat, cache.fraction = 0.1)
}

{
    NR <- 150
    NC <- 100
    mat <- RegularChunkedSparseMatrix(matrix(rpois(NR * NC, lambda=1), ncol=NC), chunks=c(5, 10))

    test_that("sparse regularly-chunked integer matrix passes basic checks", {
        expect_true(is_sparse(mat))
        expect_s4_class(chunkGrid(mat), "RegularArrayGrid")

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_true(raticate.tests::is_sparse(parsed))
        expect_true(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, cache.fraction = 0.01)
    big_test_suite(mat, cache.fraction = 0.1)
}
