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
    name <- "dense regular-chunked double matrix"

    test_that(paste(name, "passes basic checks"), {
        expect_s4_class(chunkGrid(mat), "RegularArrayGrid")
        expect_identical(type(mat), "double")

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_false(raticate.tests::sparse(parsed))
        expect_false(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, name)
}

# Repeating with a different shape, chunk size and type, just for some more thorough coverage.
{
    NR <- 75
    NC <- 50
    mat <- RegularChunkedMatrix(matrix(rpois(NR * NC, lambda=2), ncol=NC), chunks=c(5, 5))
    name <- "dense regular-chunked integer matrix"

    test_that(paste(name, "passes basic checks"), {
        expect_s4_class(chunkGrid(mat), "RegularArrayGrid")
        expect_identical(type(mat), "integer")

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_false(raticate.tests::sparse(parsed))
        expect_true(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, name)
}

{
    NR <- 0
    NC <- 10
    mat <- RegularChunkedMatrix(matrix(integer(0), nrow=NR, ncol=NC), chunks=c(NR, NC))
    name <- "dense regular-chunked integer matrix with no rows"

    test_that(paste0(name, "passes basic checks"), {
        expect_s4_class(chunkGrid(mat), "RegularArrayGrid")
        expect_identical(type(mat), "integer")

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_false(raticate.tests::sparse(parsed))
    })

    big_test_suite(mat, name)
}

{
    NR <- 10
    NC <- 0
    mat <- RegularChunkedMatrix(matrix(integer(0), nrow=NR, ncol=NC), chunks=c(NR, NC))
    name <- "dense regular-chunked integer matrix with no columns"

    test_that(paste0(name, "passes basic checks"), {
        expect_s4_class(chunkGrid(mat), "RegularArrayGrid")
        expect_identical(type(mat), "integer")

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_false(raticate.tests::sparse(parsed))
    })

    big_test_suite(mat, name)
}
