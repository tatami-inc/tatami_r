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
    paste0(prefix, "[", paste(vapply(colnames(params), function(x) paste0(x, "=", params[1,x]), ""), collapse=", "), "]")
}

simple_test_suite <- function(mat, cache.fraction) {
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
        test_that(pretty_name("dense full simple ", scenarios[i,]), {
            cache <- scenarios[i,"cache"]
            row <- scenarios[i,"row"]
            oracle <- scenarios[i,"oracle"]
            step <- scenarios[i,"step"]

            end <- if (row) nrow(mat) else ncol(mat) 
            iseq <- seq(1, end, by=step)
            if (oracle) {
                extracted <- raticate.tests::oracular_dense_full(ptr, row, iseq)
            } else {
                extracted <- raticate.tests::myopic_dense_full(ptr, row, iseq)
            }

            for (i in seq_along(iseq)) {
                j <- iseq[i]
                expected <- if (row) mat[j,] else mat[,j]
                expect_identical(extracted[[i]], expected)
            }
        })
    }
}

big_test_suite <- function(mat, cache.fraction) {
    simple_test_suite(mat, cache.fraction)
}
