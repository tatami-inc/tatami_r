library(DelayedArray)

dummy_sparse <- function(v, offset = 1L) {
    list(index = seq_along(v) + as.integer(offset) - 1L, value = v)
}

extract_sparse <- function(v, offset = 1L) {
    list(index = which(v != 0) + as.integer(offset) - 1L, value = v[v!=0])
}

get_cache_size <- function(mat, cache.fraction, sparse) {
    if (sparse) {
        # For testing, the cache type is always double + int for the indices.
        type.size <- 12
    } else {
        # For testing, the cache type is always double.
        type.size <- 8
    }
    cache.fraction * nrow(mat) * ncol(mat) * type.size
}

create_predictions <- function(iterdim, step, mode) {
    if (iterdim == 0) {
        integer(0)
    } else {
        iseq <- seq(1, iterdim, by=step)
        if (mode == "reverse") {
            rev(iseq)
        } else if (mode == "random") {
            sample(iseq)
        } else {
            iseq
        }
    }
}

pretty_name <- function(prefix, params) {
    paste0(prefix, "[", paste(vapply(colnames(params), function(x) paste0(x, "=", deparse(params[,x][[1]])), ""), collapse=", "), "]")
}

create_expected_dense <- function(mat, row, iseq, keep) {
    all.expected <- vector("list", length(iseq))
    for (i in seq_along(iseq)) {
        j <- iseq[i]
        expected <- if (row) mat[j,] else mat[,j]
        if (!is.null(keep)) {
            expected <- expected[keep]
        }
        all.expected[[i]] <- as.double(expected)
    }
    all.expected
}

fill_sparse <- function(observed, otherdim, keep) {
    for (i in seq_along(observed)) {
        both <- observed[[i]]
        if (!is.null(keep)) {
            vec <- numeric(length(keep))
            m <- match(both$index, keep)
            stopifnot(all(!is.na(m)))
            vec[m] <- both$value
        } else {
            vec <- numeric(otherdim)
            vec[both$index] <- both$value
        }
        observed[[i]] <- vec
    }
    observed
}

unlist_to_integer <- function(x) {
    if (length(x) == 0) {
        integer(0) # handle zero-length lists being unlisted.
    } else {
        unlist(x)
    }
}

full_test_suite <- function(mat) {
    scenarios <- expand.grid(
        cache = c(0, 0.01, 0.1, 0.5),
        row = c(TRUE, FALSE),
        oracle = c(FALSE, TRUE),
        mode = c("forward", "reverse", "random"), 
        step = c(1, 5),
        stringsAsFactors=FALSE
    )

    for (i in seq_len(nrow(scenarios))) {
        cache <- scenarios[i,"cache"]
        row <- scenarios[i,"row"]
        oracle <- scenarios[i,"oracle"]
        mode <- scenarios[i, "mode"]
        step <- scenarios[i,"step"]

        iterdim <- if (row) nrow(mat) else ncol(mat) 
        otherdim <- if (row) ncol(mat) else nrow(mat)
        iseq <- create_predictions(iterdim, step, mode)
        all.expected <- create_expected_dense(mat, row, iseq, NULL)

        test_that(pretty_name("dense full ", scenarios[i,]), {
            cache.size <- get_cache_size(mat, cache, sparse=FALSE)
            ptr <- raticate.tests::parse(mat, cache.size, cache.size > 0)

            if (oracle) {
                extracted <- raticate.tests::oracular_dense_full(ptr, row, iseq)
            } else {
                extracted <- raticate.tests::myopic_dense_full(ptr, row, iseq)
            }

            expect_identical(extracted, all.expected)
        })

        test_that(pretty_name("sparse full ", scenarios[i,]), {
            cache.size <- get_cache_size(mat, cache, sparse=TRUE)
            ptr <- raticate.tests::parse(mat, cache.size, cache.size > 0)

            if (oracle) {
                FUN <- raticate.tests::oracular_sparse_full
            } else {
                FUN <- raticate.tests::myopic_sparse_full
            }

            extractor.b <- FUN(ptr, row, iseq, TRUE, TRUE)
            expect_identical(all.expected, fill_sparse(extractor.b, otherdim, NULL))

            extractor.i <- FUN(ptr, row, iseq, FALSE, TRUE)
            expect_identical(extractor.i, lapply(extractor.b, function(y) y$index))

            extractor.v <- FUN(ptr, row, iseq, TRUE, FALSE)
            expect_identical(extractor.v, lapply(extractor.b, function(y) y$value))

            extractor.n <- unlist_to_integer(FUN(ptr, row, iseq, FALSE, FALSE))
            expect_identical(extractor.n, lengths(extractor.v))

            if (DelayedArray::is_sparse(mat)) {
                prod <- length(iseq) * otherdim
                if (prod > 0) {
                    expect_true(sum(extractor.n) < prod)
                }
            }
        })
    }
}

block_test_suite <- function(mat) {
    scenarios <- expand.grid(
        cache = c(0, 0.01, 0.1, 0.5),
        row = c(TRUE, FALSE),
        oracle = c(FALSE, TRUE),
        mode = c("forward", "reverse", "random"), 
        step = c(1, 5),
        block = list(c(0, 0.3), c(0.2, 0.66), c(0.6, 0.37)),
        stringsAsFactors=FALSE
    )

    for (i in seq_len(nrow(scenarios))) {
        cache <- scenarios[i,"cache"]
        row <- scenarios[i,"row"]
        oracle <- scenarios[i,"oracle"]
        mode <- scenarios[i, "mode"]
        step <- scenarios[i,"step"]
        block <- scenarios[i,"block"][[1]]

        iterdim <- if (row) nrow(mat) else ncol(mat) 
        otherdim <- if (row) ncol(mat) else nrow(mat)
        iseq <- create_predictions(iterdim, step, mode)

        bstart <- floor(block[[1]] * otherdim) + 1L
        blen <- floor(block[[2]] * otherdim)
        keep <- (bstart - 1L) + seq_len(blen)
        all.expected <- create_expected_dense(mat, row, iseq, keep)

        test_that(pretty_name("dense block ", scenarios[i,]), {
            cache.size <- get_cache_size(mat, cache, sparse=FALSE)
            ptr <- raticate.tests::parse(mat, cache.size, cache.size > 0)

            if (oracle) {
                extracted <- raticate.tests::oracular_dense_block(ptr, row, iseq, bstart, blen) 
            } else {
                extracted <- raticate.tests::myopic_dense_block(ptr, row, iseq, bstart, blen)
            }

            expect_identical(extracted, all.expected)
        })

        test_that(pretty_name("sparse block ", scenarios[i,]), {
            cache.size <- get_cache_size(mat, cache, sparse=TRUE)
            ptr <- raticate.tests::parse(mat, cache.size, cache.size > 0)

            if (oracle) {
                FUN <- raticate.tests::oracular_sparse_block
            } else {
                FUN <- raticate.tests::myopic_sparse_block
            }

            extractor.b <- FUN(ptr, row, iseq, bstart, blen, TRUE, TRUE)
            expect_identical(all.expected, fill_sparse(extractor.b, otherdim, keep))

            extractor.i <- FUN(ptr, row, iseq, bstart, blen, FALSE, TRUE)
            expect_identical(extractor.i, lapply(extractor.b, function(y) y$index))

            extractor.v <- FUN(ptr, row, iseq, bstart, blen, TRUE, FALSE)
            expect_identical(extractor.v, lapply(extractor.b, function(y) y$value))

            extractor.n <- unlist_to_integer(FUN(ptr, row, iseq, bstart, blen, FALSE, FALSE))
            expect_identical(extractor.n, lengths(extractor.v))

            if (DelayedArray::is_sparse(mat)) {
                prod <- length(iseq) * blen
                if (prod > 0) {
                    expect_true(sum(extractor.n) < prod)
                }
            }
        })
    }
}

index_test_suite <- function(mat) {
    scenarios <- expand.grid(
        cache = c(0, 0.01, 0.1, 0.5),
        row = c(TRUE, FALSE),
        oracle = c(FALSE, TRUE),
        mode = c("forward", "reverse", "random"), 
        step = c(1, 5),
        index = list(c(0, 3), c(0.33, 4), c(0.5, 5)),
        stringsAsFactors=FALSE
    )

    for (i in seq_len(nrow(scenarios))) {
        cache <- scenarios[i,"cache"]
        row <- scenarios[i,"row"]
        oracle <- scenarios[i,"oracle"]
        mode <- scenarios[i, "mode"]
        step <- scenarios[i,"step"]
        index_params <- scenarios[i,"index"][[1]]

        iterdim <- if (row) nrow(mat) else ncol(mat) 
        otherdim <- if (row) ncol(mat) else nrow(mat)
        iseq <- create_predictions(iterdim, step, mode)

        istart <- floor(index_params[[1]] * otherdim) + 1L
        if (otherdim == 0) {
            keep <- integer(0)
        } else {
            keep <- seq(istart, otherdim, by=index_params[[2]])
        }
        all.expected <- create_expected_dense(mat, row, iseq, keep)

        test_that(pretty_name("dense index ", scenarios[i,]), {
            cache.size <- get_cache_size(mat, cache, sparse=FALSE)
            ptr <- raticate.tests::parse(mat, cache.size, cache.size > 0)

            if (oracle) {
                extracted <- raticate.tests::oracular_dense_indexed(ptr, row, iseq, keep) 
            } else {
                extracted <- raticate.tests::myopic_dense_indexed(ptr, row, iseq, keep)
            }

            expect_identical(all.expected, extracted)
        })

        test_that(pretty_name("sparse index ", scenarios[i,]), {
            cache.size <- get_cache_size(mat, cache, sparse=TRUE)
            ptr <- raticate.tests::parse(mat, cache.size, cache.size > 0)

            if (oracle) {
                FUN <- raticate.tests::oracular_sparse_indexed
            } else {
                FUN <- raticate.tests::myopic_sparse_indexed
            }

            extractor.b <- FUN(ptr, row, iseq, keep, TRUE, TRUE)
            expect_identical(all.expected, fill_sparse(extractor.b, otherdim, keep))

            extractor.i <- FUN(ptr, row, iseq, keep, FALSE, TRUE)
            expect_identical(extractor.i, lapply(extractor.b, function(y) y$index))

            extractor.v <- FUN(ptr, row, iseq, keep, TRUE, FALSE)
            expect_identical(extractor.v, lapply(extractor.b, function(y) y$value))

            extractor.n <- unlist_to_integer(FUN(ptr, row, iseq, keep, FALSE, FALSE))
            expect_identical(extractor.n, lengths(extractor.v))

            if (DelayedArray::is_sparse(mat)) {
                prod <- length(iseq) * length(keep)
                if (prod > 0) {
                    expect_true(sum(extractor.n) < prod)
                }
            }
        })
    }
}

reuse_test_suite <- function(mat) {
    scenarios <- expand.grid(
        cache = c(0, 0.01, 0.1, 0.5),
        row = c(TRUE, FALSE),
        oracle = c(FALSE, TRUE),
        step = c(1, 5),
        mode = c("forward", "alternating"),
        stringsAsFactors=FALSE
    )

    for (i in seq_len(nrow(scenarios))) {
        cache <- scenarios[i,"cache"]
        row <- scenarios[i,"row"]
        oracle <- scenarios[i,"oracle"]
        step <- scenarios[i,"step"]
        mode <- scenarios[i,"mode"]

        iterdim <- if (row) nrow(mat) else ncol(mat) 
        otherdim <- if (row) ncol(mat) else nrow(mat)

        # Creating a vector of predictions where we constantly double back to
        # re-use previous elements.
        iseq <- (function() {
            predictions <- list()
            i <- 0L
            while (i < iterdim) {
                current <- i + seq_len(step * 2)
                current <- current[current <= iterdim]
                if (mode == "alternating") {
                    if (length(predictions) %% 2 == 1L) {
                        current <- rev(current)
                    }
                }
                predictions <- append(predictions, list(current))
                i <- i + step
            }
            unlist_to_integer(predictions)
        })()
        all.expected <- create_expected_dense(mat, row, iseq, NULL)

        test_that(pretty_name("dense full re-used ", scenarios[i,]), {
            cache.size <- get_cache_size(mat, cache, sparse=FALSE)
            ptr <- raticate.tests::parse(mat, cache.size, cache.size > 0)

            if (oracle) {
                extracted <- raticate.tests::oracular_dense_full(ptr, row, iseq)
            } else {
                extracted <- raticate.tests::myopic_dense_full(ptr, row, iseq)
            }

            expect_identical(all.expected, extracted)
        })

        test_that(pretty_name("sparse full re-used ", scenarios[i,]), {
            cache.size <- get_cache_size(mat, cache, sparse=TRUE)
            ptr <- raticate.tests::parse(mat, cache.size, cache.size > 0)

            if (oracle) {
                FUN <- raticate.tests::oracular_sparse_full
            } else {
                FUN <- raticate.tests::myopic_sparse_full
            }

            extractor.b <- FUN(ptr, row, iseq, TRUE, TRUE)
            expect_identical(all.expected, fill_sparse(extractor.b, otherdim, NULL))
        })
    }
}

parallel_test_suite <- function(mat) {
    for (cache in c(0, 0.01, 0.1, 0.5)) {
        refr <- Matrix::rowSums(mat)
        refc <- Matrix::colSums(mat)

        test_that("dense sums", {
            cache.size <- get_cache_size(mat, cache, sparse=FALSE)
            ptr <- raticate.tests::parse(mat, cache.size, cache.size > 0)

            expect_equal(refr, raticate.tests::myopic_dense_sums(ptr, TRUE, 1))
            expect_equal(refr, raticate.tests::oracular_dense_sums(ptr, TRUE, 1))
            expect_equal(refr, raticate.tests::myopic_dense_sums(ptr, TRUE, 3))
            expect_equal(refr, raticate.tests::oracular_dense_sums(ptr, TRUE, 3))

            expect_equal(refc, raticate.tests::myopic_dense_sums(ptr, FALSE, 1))
            expect_equal(refc, raticate.tests::oracular_dense_sums(ptr, FALSE, 1))
            expect_equal(refc, raticate.tests::myopic_dense_sums(ptr, FALSE, 3))
            expect_equal(refc, raticate.tests::oracular_dense_sums(ptr, FALSE, 3))
        })

        test_that("sparse sums", {
            cache.size <- get_cache_size(mat, cache, sparse=TRUE)
            ptr <- raticate.tests::parse(mat, cache.size, cache.size > 0)

            expect_equal(refr, raticate.tests::myopic_sparse_sums(ptr, TRUE, 1))
            expect_equal(refr, raticate.tests::oracular_sparse_sums(ptr, TRUE, 1))
            expect_equal(refr, raticate.tests::myopic_sparse_sums(ptr, TRUE, 3))
            expect_equal(refr, raticate.tests::oracular_sparse_sums(ptr, TRUE, 3))

            expect_equal(refc, raticate.tests::myopic_sparse_sums(ptr, FALSE, 1))
            expect_equal(refc, raticate.tests::oracular_sparse_sums(ptr, FALSE, 1))
            expect_equal(refc, raticate.tests::myopic_sparse_sums(ptr, FALSE, 3))
            expect_equal(refc, raticate.tests::oracular_sparse_sums(ptr, FALSE, 3))
        })
    }
}

big_test_suite <- function(mat) {
    full_test_suite(mat)
    gc(full=TRUE)

    block_test_suite(mat)
    gc(full=TRUE)

    index_test_suite(mat)
    gc(full=TRUE)

    reuse_test_suite(mat)
    gc(full=TRUE)

    parallel_test_suite(mat)
    gc(full=TRUE)
}
