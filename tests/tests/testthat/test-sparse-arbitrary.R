# This tests dense matrix extraction with arbitrary grids.
# library(testthat); source("setup.R"); source("test-sparse-arbitrary.R")

setClass("ArbitraryChunkedSparseMatrix", contains="SVT_SparseMatrix", slots=c(rowticks="integer", colticks="integer"))
setMethod("chunkGrid", "ArbitraryChunkedSparseMatrix", function(x) ArbitraryArrayGrid(list(x@rowticks, x@colticks)))
ArbitraryChunkedSparseMatrix <- function(mat, numticks) {
    rt <- sort(union(sample(nrow(mat), numticks[1]), nrow(mat)))
    ct <- sort(union(sample(ncol(mat), numticks[2]), ncol(mat)))
    spmat <- as(mat, "SVT_SparseMatrix")
    new("ArbitraryChunkedSparseMatrix", spmat, rowticks=rt, colticks=ct)
}

set.seed(200000)

{
    NR <- 54
    NC <- 201
    mat <- ArbitraryChunkedSparseMatrix(Matrix::rsparsematrix(NR, NC, 0.2), numticks=c(17, 14))

    test_that("sparse arbitrarily-chunked double matrix passes basic checks", {
        expect_s4_class(chunkGrid(mat), "ArbitraryArrayGrid")
        expect_true(is_sparse(mat))

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_true(raticate.tests::is_sparse(parsed))
        expect_true(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, cache.fraction = 0.01)
    big_test_suite(mat, cache.fraction = 0.1)
}

{
    NR <- 178
    NC <- 85 
    mat <- ArbitraryChunkedSparseMatrix(matrix(rpois(NR * NC, lambda=0.1), ncol=NC), numticks=c(20, 20))

    test_that("sparse arbitrarily-chunked integer matrix passes basic checks", {
        expect_s4_class(chunkGrid(mat), "ArbitraryArrayGrid")
        expect_true(is_sparse(mat))

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_true(raticate.tests::is_sparse(parsed))
        expect_true(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, cache.fraction = 0.01)
    big_test_suite(mat, cache.fraction = 0.1)
}
