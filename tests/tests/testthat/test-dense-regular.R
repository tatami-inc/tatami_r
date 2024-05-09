# This tests the dense matrix extraction with regular grids.
# library(testthat); source("setup.R"); source("test-dense-regular.R")

setClass("RegularChunkedMatrix", contains="matrix", slots=c(chunks="integer"))
setMethod("chunkdim", "RegularChunkedMatrix", function(x) x@chunks)
RegularChunkedMatrix <- function(mat, chunks) {
    new("RegularChunkedMatrix", mat, chunks=as.integer(chunks))
}

set.seed(150000)
{
    NR <- 23
    NC <- 104
    mat <- RegularChunkedMatrix(matrix(runif(NR * NC), ncol=NC), chunks=c(6, 4)) 

    test_that("dense regular-chunked double matrix passes basic checks", {
        expect_s4_class(chunkGrid(mat), "RegularArrayGrid")
        expect_identical(type(mat), "double")

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_false(raticate.tests::sparse(parsed))
        expect_false(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, cache.fraction = 0.01)
    big_test_suite(mat, cache.fraction = 0.1)
}

{
    NR <- 75
    NC <- 50
    mat <- RegularChunkedMatrix(matrix(rpois(NR * NC, lambda=2), ncol=NC), chunks=c(5, 5))

    test_that("dense regular-chunked integer matrix passes basic checks", {
        expect_s4_class(chunkGrid(mat), "RegularArrayGrid")
        expect_identical(type(mat), "integer")

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_false(raticate.tests::sparse(parsed))
        expect_true(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, cache.fraction = 0.01)
    big_test_suite(mat, cache.fraction = 0.1)
}
