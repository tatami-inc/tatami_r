#ifndef TATAMI_R_SIMPLEMATRIX_HPP
#define TATAMI_R_SIMPLEMATRIX_HPP

#include "tatami/tatami.hpp"
#include <algorithm>

namespace tatami_r { 

template<class InputObject_, typename Value_, typename CachedValue_>
void parse_simple_matrix_internal(const InputObject_& y, std::vector<CachedValue_>& cache, bool transpose) {
    cache.resize(y.size());
    if (transpose) {
        // y is a column-major matrix, but transpose() expects a row-major
        // input, so we just conceptually transpose it.
        tatami::transpose(static_cast<const Value_*>(y.begin()), cache.data(), y.cols(), y.rows());
    } else {
        std::copy(y.begin(), y.end(), cache.data());
    }
}

template<typename CachedValue_>
void parse_simple_matrix(const Rcpp::RObject& seed, std::vector<CachedValue_>& cache, bool transpose) {
    if (seed.sexp_type() == REALSXP) {
        Rcpp::NumericMatrix y(seed);
        parse_simple_matrix_internal<double>(y, cache, transpose);
    } else if (seed.sexp_type() == INTSXP) {
        Rcpp::IntegerMatrix y(seed);
        parse_simple_matrix_internal<int>(y, cache, transpose);
    } else if (seed.sexp_type() == LGLSXP) {
        Rcpp::LogicalMatrix y(seed);
        parse_simple_matrix_internal<int>(y, cache, transpose);
    }
}

}

#endif
