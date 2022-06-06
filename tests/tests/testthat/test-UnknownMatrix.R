# Check for correct parsing with SparseArraySeed objects.
# library(testthat); library(raticate.tests); source("test-UnknownMatrix.R")

library(Matrix)
library(DelayedArray)

dummy_sparse <- function(v, offset = 1L) {
    list(index = seq_along(v) + as.integer(offset) - 1L, value = v)
}

extract_sparse <- function(v, offset = 1L) {
    list(index = which(v != 0) + as.integer(offset) - 1L, value = v[v!=0])
}

test_that("Works for dense UnknownMatrix objects", {
    y <- DelayedArray(matrix(runif(100), 20, 5))
    y <- log1p(y) + rnorm(base::nrow(y))

    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), 20L)
    expect_identical(raticate.tests::ncol(z), 5L)

    expect_identical(raticate.tests::row(z, 11), y[11,])
    expect_identical(raticate.tests::column(z, 4), y[,4])
    expect_identical(raticate.tests::row_subset(z, 10, 2, 4), y[10,2:4])
    expect_identical(raticate.tests::column_subset(z, 3, 11, 15), y[11:15,3])

    expect_identical(raticate.tests::sparse_row(z, 9), dummy_sparse(y[9,]))
    expect_identical(raticate.tests::sparse_column(z, 3), dummy_sparse(y[,3]))
    expect_identical(raticate.tests::sparse_row_subset(z, 3, 1, 3), dummy_sparse(y[3,1:3], offset = 1))
    expect_identical(raticate.tests::sparse_column_subset(z, 5, 2, 19), dummy_sparse(y[2:19,5], offset = 2))

    expect_identical(raticate.tests::rows(z), lapply(seq_len(base::nrow(y)), function(i) y[i,]))
    expect_identical(raticate.tests::rows_subset(z, 3, 5), lapply(seq_len(base::nrow(y)), function(i) y[i,3:5]))
    expect_identical(raticate.tests::columns(z), lapply(seq_len(base::ncol(y)), function(i) y[,i]))
    expect_identical(raticate.tests::columns_subset(z, 5, 16), lapply(seq_len(base::ncol(y)), function(i) y[5:16,i]))

    expect_identical(raticate.tests::sparse_rows(z), lapply(seq_len(base::nrow(y)), function(i) dummy_sparse(y[i,])))
    expect_identical(raticate.tests::sparse_rows_subset(z, 3, 5), lapply(seq_len(base::nrow(y)), function(i) dummy_sparse(y[i,3:5], offset = 3)))
    expect_identical(raticate.tests::sparse_columns(z), lapply(seq_len(base::ncol(y)), function(i) dummy_sparse(y[,i])))
    expect_identical(raticate.tests::sparse_columns_subset(z, 5, 16), lapply(seq_len(base::ncol(y)), function(i) dummy_sparse(y[5:16,i], offset = 5)))

    # check that parallelization works correctly.
    rs <- Matrix::rowSums(y)
    expect_identical(raticate.tests::rowsums(z), rs)
    expect_identical(raticate.tests::rowsums_manual(z), rs)
})

test_that("Works for sparse UnknownMatrix objects", {
    y <- DelayedArray(abs(rsparsematrix(10, 20, 0.2)))
    y <- log1p(y) * runif(base::nrow(y))

    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), 10L)
    expect_identical(raticate.tests::ncol(z), 20L)

    expect_identical(raticate.tests::row(z, 6), y[6,])
    expect_identical(raticate.tests::column(z, 13), y[,13])
    expect_identical(raticate.tests::row_subset(z, 10, 3, 15), y[10,3:15])
    expect_identical(raticate.tests::column_subset(z, 18, 3, 8), y[3:8,18])

    expect_identical(raticate.tests::sparse_row(z, 7), extract_sparse(y[7,]))
    expect_identical(raticate.tests::sparse_column(z, 14), extract_sparse(y[,14]))
    expect_identical(raticate.tests::sparse_row_subset(z, 5, 2, 8), extract_sparse(y[5,2:8], offset = 2))
    expect_identical(raticate.tests::sparse_column_subset(z, 6, 3, 10), extract_sparse(y[3:10,6], offset = 3))

    expect_identical(raticate.tests::rows(z), lapply(seq_len(base::nrow(y)), function(i) y[i,]))
    expect_identical(raticate.tests::rows_subset(z, 2, 17), lapply(seq_len(base::nrow(y)), function(i) y[i,2:17]))
    expect_identical(raticate.tests::columns(z), lapply(seq_len(base::ncol(y)), function(i) y[,i]))
    expect_identical(raticate.tests::columns_subset(z, 3, 8), lapply(seq_len(base::ncol(y)), function(i) y[3:8,i]))

    expect_identical(raticate.tests::sparse_rows(z), lapply(seq_len(base::nrow(y)), function(i) extract_sparse(y[i,])))
    expect_identical(raticate.tests::sparse_rows_subset(z, 3, 15), lapply(seq_len(base::nrow(y)), function(i) extract_sparse(y[i,3:15], offset = 3)))
    expect_identical(raticate.tests::sparse_columns(z), lapply(seq_len(base::ncol(y)), function(i) extract_sparse(y[,i])))
    expect_identical(raticate.tests::sparse_columns_subset(z, 5, 7), lapply(seq_len(base::ncol(y)), function(i) extract_sparse(y[5:7,i], offset = 5)))

    rs <- Matrix::rowSums(y)
    expect_identical(raticate.tests::rowsums(z), rs)
    expect_identical(raticate.tests::rowsums_manual(z), rs)
})

test_that("Behaves correctly with small block sizes", {
    y <- DelayedArray(abs(rsparsematrix(103, 51, 0.2))) # slightly offset from multiple, check for correct capping by dimension.
    y <- log1p(y) / rnorm(base::nrow(y))

    blocksize <- 500
    DelayedArray::setAutoBlockSize(blocksize * 8)
    expect_equal(base::ncol(colAutoGrid(y)), ceiling(base::ncol(y) / floor(blocksize/base::nrow(y))))
    expect_equal(base::nrow(rowAutoGrid(y)), ceiling(base::nrow(y) / floor(blocksize/base::ncol(y))))

    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), 103L)
    expect_identical(raticate.tests::ncol(z), 51L)

    expect_identical(raticate.tests::rows(z), lapply(seq_len(base::nrow(y)), function(i) y[i,]))
    expect_identical(raticate.tests::rows_subset(z, 3, 5), lapply(seq_len(base::nrow(y)), function(i) y[i,3:5]))
    expect_identical(raticate.tests::columns(z), lapply(seq_len(base::ncol(y)), function(i) y[,i]))
    expect_identical(raticate.tests::columns_subset(z, 5, 16), lapply(seq_len(base::ncol(y)), function(i) y[5:16,i]))

    expect_identical(raticate.tests::sparse_rows(z), lapply(seq_len(base::nrow(y)), function(i) extract_sparse(y[i,])))
    expect_identical(raticate.tests::sparse_rows_subset(z, 3, 15), lapply(seq_len(base::nrow(y)), function(i) extract_sparse(y[i,3:15], offset = 3)))
    expect_identical(raticate.tests::sparse_columns(z), lapply(seq_len(base::ncol(y)), function(i) extract_sparse(y[,i])))
    expect_identical(raticate.tests::sparse_columns_subset(z, 5, 7), lapply(seq_len(base::ncol(y)), function(i) extract_sparse(y[5:7,i], offset = 5)))

    rs <- Matrix::rowSums(y)
    expect_identical(raticate.tests::rowsums(z), rs)
    expect_identical(raticate.tests::rowsums_manual(z), rs)

    setAutoBlockSize()
})

test_that("Behaves correctly with small chunk sizes", {
    setClass("MyOwnMatrix", contains="dgCMatrix")
    setMethod("chunkdim", "MyOwnMatrix", function(x) c(11L, 9L))

    y <- DelayedArray(new("MyOwnMatrix", rsparsematrix(39, 55, 0.2))) 
    y <- y * 2

    # Induce some interesting interaction with the block size.
    blocksize <- 500
    DelayedArray::setAutoBlockSize(blocksize * 8)
    expect_true(base::ncol(colAutoGrid(y)) > 1)
    expect_true(base::nrow(rowAutoGrid(y)) > 1)

    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), 39L)
    expect_identical(raticate.tests::ncol(z), 55L)

    expect_identical(raticate.tests::rows(z), lapply(seq_len(base::nrow(y)), function(i) y[i,]))
    expect_identical(raticate.tests::rows_subset(z, 3, 5), lapply(seq_len(base::nrow(y)), function(i) y[i,3:5]))
    expect_identical(raticate.tests::columns(z), lapply(seq_len(base::ncol(y)), function(i) y[,i]))
    expect_identical(raticate.tests::columns_subset(z, 5, 16), lapply(seq_len(base::ncol(y)), function(i) y[5:16,i]))

    expect_identical(raticate.tests::sparse_rows(z), lapply(seq_len(base::nrow(y)), function(i) extract_sparse(y[i,])))
    expect_identical(raticate.tests::sparse_rows_subset(z, 3, 15), lapply(seq_len(base::nrow(y)), function(i) extract_sparse(y[i,3:15], offset = 3)))
    expect_identical(raticate.tests::sparse_columns(z), lapply(seq_len(base::ncol(y)), function(i) extract_sparse(y[,i])))
    expect_identical(raticate.tests::sparse_columns_subset(z, 5, 7), lapply(seq_len(base::ncol(y)), function(i) extract_sparse(y[5:7,i], offset = 5)))

    rs <- Matrix::rowSums(y)
    expect_identical(raticate.tests::rowsums(z), rs)
    expect_identical(raticate.tests::rowsums_manual(z), rs)

    setAutoBlockSize()
})
