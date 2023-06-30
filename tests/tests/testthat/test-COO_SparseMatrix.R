# Check for correct parsing with COO_SparseMatrix objects.
# library(testthat); library(raticate.tests); source("setup.R"); source("test-COO_SparseMatrix.R")

library(Matrix)
library(DelayedArray)

test_that("Works for numeric COO_SparseMatrix objects", {
    x <- rsparsematrix(100, 20, 0.05)
    y <- as(x, "COO_SparseMatrix")

    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), 100L)
    expect_identical(raticate.tests::ncol(z), 20L)

    expect_identical(raticate.tests::row(z, 5), x[5,])
    expect_identical(raticate.tests::column(z, 9), x[,9])
})

test_that("Works for integer COO_SparseMatrix objects", {
    x <- rsparsematrix(40, 50, 0.05)
    y <- as(x, "COO_SparseMatrix")
    y@nzdata <- as.integer(ceiling(y@nzdata * 10))

    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), 40L)
    expect_identical(raticate.tests::ncol(z), 50L)

    ref <- DelayedArray(y)
    expect_identical(raticate.tests::row(z, 7), as.double(ref[7,]))
    expect_identical(raticate.tests::column(z, 2), as.double(ref[,2]))
})

test_that("Works for logical COO_SparseMatrix objects", {
    x <- rsparsematrix(25, 25, 0.05) != 0
    y <- as(x, "COO_SparseMatrix")

    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), 25L)
    expect_identical(raticate.tests::ncol(z), 25L)

    expect_identical(raticate.tests::row(z, 10), as.double(x[10,]))
    expect_identical(raticate.tests::column(z, 10), as.double(x[,10]))
})

test_that("Handles CSC COO_SparseMatrix objects in sparse mode", {
    x <- rsparsematrix(100, 50, 0.05)
    y <- as(x, "COO_SparseMatrix")
    o <- order(y@nzcoo[,2], y@nzcoo[,1])
    y@nzcoo <- y@nzcoo[o,]
    y@nzdata <- y@nzdata[o]

    # CSR'ness is preserved.
    out <- extract_sparse_array(y, list(NULL, NULL))
    expect_s4_class(out, "COO_SparseMatrix")
    expect_identical(nzcoo(out), nzcoo(out))

    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), base::nrow(x))
    expect_identical(raticate.tests::ncol(z), base::ncol(x))

    for (i in seq_len(base::nrow(x))) {
        expect_identical(raticate.tests::sparse_row(z, i), extract_sparse(x[i,]))
    }

    for (i in seq_len(base::ncol(x))) {
        expect_identical(raticate.tests::sparse_column(z, i), extract_sparse(x[,i]))
    }
})

test_that("Handles CSR COO_SparseMatrix objects", {
    x <- rsparsematrix(90, 50, 0.2)
    y <- as(x, "COO_SparseMatrix")
    o <- order(y@nzcoo[,1], y@nzcoo[,2])
    y@nzcoo <- y@nzcoo[o,]
    y@nzdata <- y@nzdata[o]

    # CSR'ness is preserved.
    out <- extract_sparse_array(y, list(NULL, NULL))
    expect_s4_class(out, "COO_SparseMatrix")
    expect_identical(nzcoo(out), nzcoo(out))

    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), base::nrow(x))
    expect_identical(raticate.tests::ncol(z), base::ncol(x))

    for (i in seq_len(base::nrow(x))) {
        expect_identical(raticate.tests::sparse_row(z, i), extract_sparse(x[i,]))
    }

    for (i in seq_len(base::ncol(x))) {
        expect_identical(raticate.tests::sparse_column(z, i), extract_sparse(x[,i]))
    }
})

test_that("Handles unsorted COO_SparseMatrix objects", {
    x <- rsparsematrix(67, 80, 0.1)
    y <- as(x, "COO_SparseMatrix")
    o <- sample(base::nrow(y@nzcoo))
    y@nzcoo <- y@nzcoo[o,]
    y@nzdata <- y@nzdata[o]

    # Unsortedness is preserved.
    out <- extract_sparse_array(y, list(NULL, NULL))
    expect_s4_class(out, "COO_SparseMatrix")
    expect_identical(nzcoo(out), nzcoo(out))

    z <- raticate.tests::parse(y)
    expect_identical(raticate.tests::nrow(z), base::nrow(x))
    expect_identical(raticate.tests::ncol(z), base::ncol(x))

    for (i in seq_len(base::nrow(x))) {
        expect_identical(raticate.tests::sparse_row(z, i), extract_sparse(x[i,]))
    }

    for (i in seq_len(base::ncol(x))) {
        expect_identical(raticate.tests::sparse_column(z, i), extract_sparse(x[,i]))
    }
})
