# library(testthat); source("setup.R"); source("test-sparse-regular.R")

setClass("RegularChunkedSparseMatrix", contains="SVT_SparseMatrix", slots=c(chunks="integer"))
setMethod("chunkdim", "RegularChunkedSparseMatrix", function(x) x@chunks)
RegularChunkedSparseMatrix <- function(mat, chunks) {
    spmat <- as(mat, "SVT_SparseMatrix")
    new("RegularChunkedSparseMatrix", spmat, chunks=as.integer(chunks))
}

set.seed(150000)

{
    NR <- 24
    NC <- 104
    mat <- RegularChunkedSparseMatrix(Matrix::rsparsematrix(NR, NC, density=0.24), chunks=c(8, 7)) 
    name <- "sparse regularly-chunked double matrix"

    test_that(paste(name, "passes basic checks"), {
        expect_true(is_sparse(mat))
        expect_s4_class(chunkGrid(mat), "RegularArrayGrid")

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_true(raticate.tests::sparse(parsed))
        expect_false(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, name)
}

{
    NR <- 75
    NC <- 50
    mat <- RegularChunkedSparseMatrix(matrix(rpois(NR * NC, lambda=1), ncol=NC), chunks=c(3, 10))
    name <- "sparse regularly-chunked integer matrix"

    test_that(paste(name, "passes basic checks"), {
        expect_true(is_sparse(mat))
        expect_s4_class(chunkGrid(mat), "RegularArrayGrid")

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_true(raticate.tests::sparse(parsed))
        expect_true(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, name)
}

{
    NR <- 0
    NC <- 10
    mat <- RegularChunkedSparseMatrix(matrix(integer(0), nrow=NR, ncol=NC), chunks=c(NR, NC))
    name <- "sparse regular-chunked integer matrix with no rows" 

    test_that(paste(name, "passes basic checks"), {
        expect_s4_class(chunkGrid(mat), "RegularArrayGrid")
        expect_identical(type(mat), "integer")

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_true(raticate.tests::sparse(parsed))
    })

    big_test_suite(mat, name)
}

{
    NR <- 10
    NC <- 0
    mat <- RegularChunkedSparseMatrix(matrix(integer(0), nrow=NR, ncol=NC), chunks=c(NR, NC))
    name <- "sparse regular-chunked integer matrix with no columns" 

    test_that(paste(name, "passes basic checks"), {
        expect_s4_class(chunkGrid(mat), "RegularArrayGrid")
        expect_identical(type(mat), "integer")

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_true(raticate.tests::sparse(parsed))
    })

    big_test_suite(mat, name)
}
