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
    name <- "dense arbitrary-chunked double matrix"

    test_that(paste(name, "passes basic checks"), {
        expect_type(mat, "double")
        expect_s4_class(chunkGrid(mat), "ArbitraryArrayGrid")

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_false(raticate.tests::sparse(parsed))
        expect_false(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, name)
}

# Repeating with a different shape, chunk size and type, just for some more thorough coverage.
{
    NR <- 97
    NC <- 46
    mat <- ArbitraryChunkedMatrix(matrix(rpois(NR * NC, lambda=10), ncol=NC), numticks=c(19, 15))
    name <- "dense arbitrary-chunked integer matrix"

    test_that(paste(name, "passes basic checks"), {
        expect_type(mat, "integer")
        expect_s4_class(chunkGrid(mat), "ArbitraryArrayGrid")

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_false(raticate.tests::sparse(parsed))
        expect_true(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, name)
}

{
    NR <- 0
    NC <- 10
    mat <- ArbitraryChunkedMatrix(matrix(double(0), nrow=NR, ncol=NC), numticks=c(NR, NC))
    name <- "dense arbitrary-chunked double matrix with no rows"

    test_that(paste(name, "passes basic checks"), {
        expect_type(mat, "double")
        expect_s4_class(chunkGrid(mat), "ArbitraryArrayGrid")

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_false(raticate.tests::sparse(parsed))
    })

    big_test_suite(mat, name)
}

{
    NR <- 10
    NC <- 0
    mat <- ArbitraryChunkedMatrix(matrix(double(0), nrow=NR, ncol=NC), numticks=c(NR, NC))
    name <- "dense arbitrary-chunked double matrix with no columns"

    test_that(paste(name, "passes basic checks"), {
        expect_type(mat, "double")
        expect_s4_class(chunkGrid(mat), "ArbitraryArrayGrid")

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_false(raticate.tests::sparse(parsed))
    })

    big_test_suite(mat, name)
}
