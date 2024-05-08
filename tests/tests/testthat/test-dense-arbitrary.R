# This tests dense matrix extraction with arbitrary grids.
# library(testthat); source("setup.R"); source("test-dense-arbitrary.R")

setClass("ArbitraryChunkedMatrix", contains="matrix", slots=c(rowticks="integer", colticks="integer"))
setMethod("chunkGrid", "ArbitraryChunkedMatrix", function(x) ArbitraryArrayGrid(list(x@rowticks, x@colticks)))
ArbitraryChunkedMatrix <- function(mat, numticks) {
    rt <- sort(union(sample(nrow(mat), numticks[1]), nrow(mat)))
    ct <- sort(union(sample(ncol(mat), numticks[2]), ncol(mat)))
    new("ArbitraryChunkedMatrix", mat, rowticks=rt, colticks=ct)
}

set.seed(200000)

{
    NR <- 31
    NC <- 89 
    mat <- ArbitraryChunkedMatrix(matrix(runif(NR * NC), ncol=NC), numticks=c(11L, 20L))

    test_that("dense arbitrary-chunked double matrix passes basic checks", {
        expect_type(mat, "double")
        expect_s4_class(chunkGrid(mat), "ArbitraryArrayGrid")

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_false(raticate.tests::is_sparse(parsed))
        expect_false(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, cache.fraction = 0.01)
    big_test_suite(mat, cache.fraction = 0.1)
}

{
    NR <- 97
    NC <- 46
    mat <- ArbitraryChunkedMatrix(matrix(rpois(NR * NC, lambda=2), ncol=NC), numticks=c(19, 15))

    test_that("dense arbitrary-chunked integer matrix passes basic checks", {
        expect_type(mat, "integer")
        expect_s4_class(chunkGrid(mat), "ArbitraryArrayGrid")

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_false(raticate.tests::is_sparse(parsed))
        expect_true(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, cache.fraction = 0.01)
    big_test_suite(mat, cache.fraction = 0.1)
}
