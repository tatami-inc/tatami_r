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
    NR <- 27
    NC <- 101
    mat <- ArbitraryChunkedSparseMatrix(Matrix::rsparsematrix(NR, NC, 0.2), numticks=c(17, 14))
    name <- "sparse arbitrarily-chunked double matrix"

    test_that(paste(name, "passes basic checks"), {
        expect_s4_class(chunkGrid(mat), "ArbitraryArrayGrid")
        expect_true(is_sparse(mat))

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_true(raticate.tests::sparse(parsed))
        expect_true(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, name)
}

{
    NR <- 89
    NC <- 43 
    mat <- ArbitraryChunkedSparseMatrix(matrix(rpois(NR * NC, lambda=0.1), ncol=NC), numticks=c(20, 21))
    name <- "sparse arbitrarily-chunked integer matrix"

    test_that(paste(name, "passes basic checks"), {
        expect_s4_class(chunkGrid(mat), "ArbitraryArrayGrid")
        expect_true(is_sparse(mat))

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_true(raticate.tests::sparse(parsed))
        expect_false(raticate.tests::prefer_rows(parsed))
    })

    big_test_suite(mat, name)
}

{
    NR <- 0
    NC <- 10
    mat <- ArbitraryChunkedSparseMatrix(matrix(double(0), nrow=NR, ncol=NC), numticks=c(NR, NC))
    name <- "sparse arbitrarily-chunked double matrix with no rows"

    test_that(paste(name, "passes basic checks"), {
        expect_s4_class(chunkGrid(mat), "ArbitraryArrayGrid")
        expect_true(is_sparse(mat))

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_true(raticate.tests::sparse(parsed))
    })

    big_test_suite(mat, name)
}

{
    NR <- 10
    NC <- 0
    mat <- ArbitraryChunkedSparseMatrix(matrix(double(0), nrow=NR, ncol=NC), numticks=c(NR, NC))
    name <- "sparse arbitrarily-chunked double matrix with no columns"

    test_that(paste(name, "passes basic checks"), {
        expect_s4_class(chunkGrid(mat), "ArbitraryArrayGrid")
        expect_true(is_sparse(mat))

        parsed <- raticate.tests::parse(mat, 0, FALSE)
        expect_true(raticate.tests::sparse(parsed))
    })

    big_test_suite(mat, name)
}
