dummy_sparse <- function(v, offset = 1L) {
    list(index = seq_along(v) + as.integer(offset) - 1L, value = v)
}

extract_sparse <- function(v, offset = 1L) {
    list(index = which(v != 0) + as.integer(offset) - 1L, value = v[v!=0])
}

get_cache_size <- function(mat, cache.fraction) {
    type.size <- c(logical = 4L, integer = 4L, numeric = 8L, double = 8L)[[DelayedArray::type(mat)]]
    cache.fraction * nrow(mat) * ncol(mat) * type.size
}

pretty_name <- function(prefix, params) {
    paste0(prefix, "[", paste(vapply(colnames(params), function(x) paste0(x, "=", deparse(params[,x][[1]])), ""), collapse=", "), "]")
}

full_test_suite <- function(mat, cache.fraction) {
    cache.size <- get_cache_size(mat, cache.fraction)
    ptr <- raticate.tests::parse(mat, cache.size, cache.size > 0)
    expect_identical(nrow(mat), raticate.tests::nrow(ptr))
    expect_identical(ncol(mat), raticate.tests::ncol(ptr))

    scenarios <- expand.grid(
        cache = cache.fraction,
        row = c(TRUE, FALSE),
        oracle = c(FALSE, TRUE),
        step = c(1, 5, 10)
    )

    for (i in seq_len(nrow(scenarios))) {
        cache <- scenarios[i,"cache"]
        row <- scenarios[i,"row"]
        oracle <- scenarios[i,"oracle"]
        step <- scenarios[i,"step"]

        iterdim <- if (row) nrow(mat) else ncol(mat) 
        otherdim <- if (row) ncol(mat) else nrow(mat)
        iseq <- seq(1, iterdim, by=step)

        test_that(pretty_name("dense full simple ", scenarios[i,]), {
            if (oracle) {
                extracted <- raticate.tests::oracular_dense_full(ptr, row, iseq)
            } else {
                extracted <- raticate.tests::myopic_dense_full(ptr, row, iseq)
            }
            for (i in seq_along(iseq)) {
                j <- iseq[i]
                expected <- if (row) mat[j,] else mat[,j]
                expected <- as.double(expected)
                expect_identical(extracted[[i]], expected)
            }
        })

        test_that(pretty_name("sparse full simple ", scenarios[i,]), {
            if (oracle) {
                FUN <- raticate.tests::oracular_sparse_full
            } else {
                FUN <- raticate.tests::myopic_sparse_full
            }
            extractor.b <- FUN(ptr, row, iseq, TRUE, TRUE)
            extractor.i <- FUN(ptr, row, iseq, FALSE, TRUE)
            extractor.v <- FUN(ptr, row, iseq, TRUE, FALSE)
            extractor.n <- FUN(ptr, row, iseq, FALSE, FALSE)

            ncount <- 0L
            for (i in seq_along(iseq)) {
                j <- iseq[i]
                expected <- if (row) mat[j,] else mat[,j]
                expected <- as.double(expected)

                both <- extractor.b[[i]]
                observed <- numeric(otherdim)
                observed[both$index] <- both$value
                expect_identical(observed, expected)

                expect_identical(both$index, extractor.i[[i]])
                expect_identical(both$value, extractor.v[[i]])
                expect_identical(length(both$value), extractor.n[[i]])

                ncount <- ncount + length(both$value)
            }

            if (DelayedArray::is_sparse(mat)) {
                expect_true(ncount < nrow(mat) * ncol(mat))
            }
        })
    }
}

block_test_suite <- function(mat, cache.fraction) {
    cache.size <- get_cache_size(mat, cache.fraction)
    ptr <- raticate.tests::parse(mat, cache.size, cache.size > 0)
    expect_identical(nrow(mat), raticate.tests::nrow(ptr))
    expect_identical(ncol(mat), raticate.tests::ncol(ptr))

    scenarios <- expand.grid(
        cache = cache.fraction,
        row = c(TRUE, FALSE),
        oracle = c(FALSE, TRUE),
        step = c(1, 5, 10),
        block = list(c(0, 0.3), c(0.2, 0.66), c(0.6, 0.37))
    )

    for (i in seq_len(nrow(scenarios))) {
        cache <- scenarios[i,"cache"]
        row <- scenarios[i,"row"]
        oracle <- scenarios[i,"oracle"]
        step <- scenarios[i,"step"]

        iterdim <- if (row) nrow(mat) else ncol(mat) 
        otherdim <- if (row) ncol(mat) else nrow(mat)
        iseq <- seq(1, iterdim, by=step)

        block <- scenarios[i,"block"][[1]]
        bstart <- floor(block[[1]] * otherdim) + 1L
        blen <- floor(block[[2]] * otherdim)
        keep <- (bstart - 1L) + seq_len(blen)

        test_that(pretty_name("dense block simple ", scenarios[i,]), {
            if (oracle) {
                extracted <- raticate.tests::oracular_dense_block(ptr, row, iseq, bstart, blen) 
            } else {
                extracted <- raticate.tests::myopic_dense_block(ptr, row, iseq, bstart, blen)
            }
            for (i in seq_along(iseq)) {
                j <- iseq[i]
                expected <- if (row) mat[j,keep] else mat[keep,j]
                expected <- as.double(expected)
                expect_identical(extracted[[i]], expected)
            }
        })

        test_that(pretty_name("sparse block simple ", scenarios[i,]), {
            if (oracle) {
                FUN <- raticate.tests::oracular_sparse_block
            } else {
                FUN <- raticate.tests::myopic_sparse_block
            }
            extractor.b <- FUN(ptr, row, iseq, bstart, blen, TRUE, TRUE)
            extractor.i <- FUN(ptr, row, iseq, bstart, blen, FALSE, TRUE)
            extractor.v <- FUN(ptr, row, iseq, bstart, blen, TRUE, FALSE)
            extractor.n <- FUN(ptr, row, iseq, bstart, blen, FALSE, FALSE)

            ncount <- 0L
            for (i in seq_along(iseq)) {
                j <- iseq[i]
                expected <- if (row) mat[j,keep] else mat[keep,j]
                expected <- as.double(expected)

                both <- extractor.b[[i]]
                observed <- numeric(otherdim)
                observed[both$index] <- both$value
                expect_identical(observed[keep], expected)

                expect_identical(both$index, extractor.i[[i]])
                expect_identical(both$value, extractor.v[[i]])
                expect_identical(length(both$value), extractor.n[[i]])

                ncount <- ncount + length(both$value)
            }

            if (DelayedArray::is_sparse(mat)) {
                expect_true(ncount < nrow(mat) * ncol(mat))
            }
        })
    }
}

index_test_suite <- function(mat, cache.fraction) {
    cache.size <- get_cache_size(mat, cache.fraction)
    ptr <- raticate.tests::parse(mat, cache.size, cache.size > 0)
    expect_identical(nrow(mat), raticate.tests::nrow(ptr))
    expect_identical(ncol(mat), raticate.tests::ncol(ptr))

    scenarios <- expand.grid(
        cache = cache.fraction,
        row = c(TRUE, FALSE),
        oracle = c(FALSE, TRUE),
        step = c(1, 5, 10),
        index = list(c(0, 3), c(0.33, 4), c(0.5, 5))
    )

    for (i in seq_len(nrow(scenarios))) {
        cache <- scenarios[i,"cache"]
        row <- scenarios[i,"row"]
        oracle <- scenarios[i,"oracle"]
        step <- scenarios[i,"step"]

        iterdim <- if (row) nrow(mat) else ncol(mat) 
        otherdim <- if (row) ncol(mat) else nrow(mat)
        iseq <- seq(1, iterdim, by=step)

        index_params <- scenarios[i,"index"][[1]]
        istart <- floor(index_params[[1]] * otherdim) + 1L
        keep <- seq(istart, otherdim, by=index_params[[2]])

        test_that(pretty_name("dense index simple ", scenarios[i,]), {
            if (oracle) {
                extracted <- raticate.tests::oracular_dense_indexed(ptr, row, iseq, keep) 
            } else {
                extracted <- raticate.tests::myopic_dense_indexed(ptr, row, iseq, keep)
            }
            for (i in seq_along(iseq)) {
                j <- iseq[i]
                expected <- if (row) mat[j,keep] else mat[keep,j]
                expected <- as.double(expected)
                expect_identical(extracted[[i]], expected)
            }
        })

        test_that(pretty_name("sparse index simple ", scenarios[i,]), {
            if (oracle) {
                FUN <- raticate.tests::oracular_sparse_indexed
            } else {
                FUN <- raticate.tests::myopic_sparse_indexed
            }
            extractor.b <- FUN(ptr, row, iseq, keep, TRUE, TRUE)
            extractor.i <- FUN(ptr, row, iseq, keep, FALSE, TRUE)
            extractor.v <- FUN(ptr, row, iseq, keep, TRUE, FALSE)
            extractor.n <- FUN(ptr, row, iseq, keep, FALSE, FALSE)

            ncount <- 0L
            for (i in seq_along(iseq)) {
                j <- iseq[i]
                expected <- if (row) mat[j,keep] else mat[keep,j]
                expected <- as.double(expected)

                both <- extractor.b[[i]]
                observed <- numeric(otherdim)
                observed[both$index] <- both$value
                expect_identical(observed[keep], expected)

                expect_identical(both$index, extractor.i[[i]])
                expect_identical(both$value, extractor.v[[i]])
                expect_identical(length(both$value), extractor.n[[i]])

                ncount <- ncount + length(both$value)
            }

            if (DelayedArray::is_sparse(mat)) {
                expect_true(ncount < nrow(mat) * ncol(mat))
            }
        })
    }
}

big_test_suite <- function(mat, cache.fraction) {
    full_test_suite(mat, cache.fraction)
    block_test_suite(mat, cache.fraction)
    index_test_suite(mat, cache.fraction)
}
