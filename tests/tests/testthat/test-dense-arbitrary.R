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
NR <- 57
NC <- 208
mat <- ArbitraryChunkedMatrix(matrix(runif(NR * NC), ncol=NC), numticks=c(11L, 20L))
expect_s4_class(chunkGrid(mat), "ArbitraryArrayGrid")
big_test_suite(mat, cache.fraction = 0.01)
big_test_suite(mat, cache.fraction = 0.1)

NR <- 150
NC <- 100
raw <- matrix(rpois(NR * NC, lambda=2), ncol=NC)
mat <- DelayedArray(new("ArbitraryChunkedMatrix", raw, chunks=c(13, 15)))
big_test_suite(mat, cache.fraction = 0.01)
big_test_suite(mat, cache.fraction = 0.1)
