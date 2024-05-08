# This tests the dense matrix extraction.
# library(testthat); source("setup.R"); source("test-dense-regular.R")

library(DelayedArray)
setClass("RegularChunkedMatrix", contains="matrix", slots=c(chunks="integer"))
setMethod("chunkdim", "RegularChunkedMatrix", function(x) x@chunks)

set.seed(100000)
NR <- 57
NC <- 208
raw <- matrix(runif(NR * NC), ncol=NC)
mat <- DelayedArray(new("RegularChunkedMatrix", raw, chunks=c(11L, 7L)))
expect_s4_class(chunkGrid(mat), "RegularArrayGrid")
big_test_suite(mat, cache.fraction = 0.01)
big_test_suite(mat, cache.fraction = 0.1)

NR <- 150
NC <- 100
raw <- matrix(rpois(NR * NC, lambda=2), ncol=NC)
mat <- DelayedArray(new("RegularChunkedMatrix", raw, chunks=c(5L, 10L)))
big_test_suite(mat, cache.fraction = 0.01)
big_test_suite(mat, cache.fraction = 0.1)
