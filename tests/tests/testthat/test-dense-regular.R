# This tests the dense matrix extraction with regular grids.
# library(testthat); source("setup.R"); source("test-dense-regular.R")

setClass("RegularChunkedMatrix", contains="matrix", slots=c(chunks="integer"))
setMethod("chunkdim", "RegularChunkedMatrix", function(x) x@chunks)
RegularChunkedMatrix <- function(mat, chunks) {
    new("RegularChunkedMatrix", mat, chunks=as.integer(chunks))
}

set.seed(150000)
NR <- 57
NC <- 208
mat <- RegularChunkedMatrix(matrix(runif(NR * NC), ncol=NC), chunks=c(11, 7)) 
expect_s4_class(chunkGrid(mat), "RegularArrayGrid")
big_test_suite(mat, cache.fraction = 0.01)
big_test_suite(mat, cache.fraction = 0.1)

NR <- 150
NC <- 100
raw <- RegularChunkedMatrix(matrix(rpois(NR * NC, lambda=2), ncol=NC), chunks=c(5, 10))
big_test_suite(mat, cache.fraction = 0.01)
big_test_suite(mat, cache.fraction = 0.1)
