dummy_sparse <- function(v, offset = 1L) {
    list(index = seq_along(v) + as.integer(offset) - 1L, value = v)
}

extract_sparse <- function(v, offset = 1L) {
    list(index = which(v != 0) + as.integer(offset) - 1L, value = v[v!=0])
}
