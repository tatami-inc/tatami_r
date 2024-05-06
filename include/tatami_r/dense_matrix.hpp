#ifndef TATAMI_R_DENSE_MATRIX_HPP
#define TATAMI_R_DENSE_MATRIX_HPP

#include "tatami/tatami.hpp"
#include <algorithm>

namespace tatami_r { 

template<bool transpose_, typename InputValue_, class InputObject_,  typename CachedValue_>
void parse_dense_matrix_internal(const InputObject_& y, std::vector<CachedValue_>& cache, size_t start_row, size_t start_col, size_t num_rows, size_t num_cols) {
    cache.resize(num_rows * num_cols);
    auto input = static_cast<const InputValue_*>(y.begin()) + start_row + start_col * static_cast<size_t>(y.rows());
    auto output = cache.data();

    if constexpr(transpose_) {
        // y is a column-major matrix, but transpose() expects a row-major
        // input, so we just conceptually transpose it.
        tatami::transpose(input, num_cols, num_rows, y.rows(), output, num_cols);
    } else {
        for (size_t c = 0; c < num_cols; ++c) {
            std::copy_n(input, num_rows, output);
            input += y.rows();
            output += num_rows;
        }
    }
}

template<bool transpose_, typename CachedValue_>
void parse_dense_matrix(const Rcpp::RObject& seed, std::vector<CachedValue_>& cache, size_t start_row, size_t start_col, size_t num_rows, size_t num_cols) {
    if (seed.sexp_type() == REALSXP) {
        Rcpp::NumericMatrix y(seed);
        parse_dense_matrix_internal<transpose_, double>(y, cache, start_row, start_col, num_rows, num_cols);
    } else if (seed.sexp_type() == INTSXP) {
        Rcpp::IntegerMatrix y(seed);
        parse_dense_matrix_internal<transpose_, int>(y, cache, start_row, start_col, num_rows, num_cols) 
    } else if (seed.sexp_type() == LGLSXP) {
        Rcpp::LogicalMatrix y(seed);
        parse_dense_matrix_internal<transpose_, int>(y, cache, start_row, start_col, num_rows, num_cols);
    }

    auto ctype = get_class_name(seed);
    throw std::runtime_error("unsupported type '" + type + "' for a " + ctype + "object");
}

}

#endif
