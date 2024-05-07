#ifndef TATAMI_R_SPARSE_MATRIX_HPP
#define TATAMI_R_SPARSE_MATRIX_HPP

#include "utils.hpp"
#include "tatami/tatami.hpp"
#include <type_traits>

namespace tatami_r { 

template<bool transpose_, class InputObject_, SEXPTYPE desired_sexp_, typename InputValue_, typename CachedValue_, typename CachedIndex_, typename Index_>
void parse_sparse_matrix_internal(
    Rcpp::RObject seed, 
    std::vector<CachedValue_*>& value_ptrs, 
    std::vector<CachedIndex_*>& index_ptrs, 
    std::vector<size_t>& counts, 
    const Rcpp::IntegerVector& secondary_extract)
{
    auto dims = parse_dims(seed.slot("dim"));
    int NR = dims.first;
    int NC = dims.second;

    Rcpp::List svt = seed.slot("SVT");
    if (svt.size() != NC) {
        auto ctype = get_class_name(seed);
        throw std::runtime_error(std::string("'SVT' slot in a ") + ctype + " object should have length equal to the number of columns");
    }

    bool needs_value = !value_ptrs.empty();
    bool needs_index = !index_ptrs.empty();

    for (int c = 0; c < NC; ++c) {
        Rcpp::RObject raw_inner(svt[c]);
        if (raw_inner == R_NilValue) {
            continue;
        }

        Rcpp::List inner(raw_inner);
        if (inner.size() != 2) {
            auto ctype = get_class_name(seed);
            throw std::runtime_error("each entry of the 'SVT' slot of a " + ctype + " object should be a list of length 2 or NULL");
        }

        // Verify type to ensure that we're not making a view on a temporary array.
        Rcpp::RObject first = inner[0];
        if (first.sexp_type() != INTSXP) {
            auto ctype = get_class_name(seed);
            throw std::runtime_error("first entry of each element of the 'SVT' slot in a " + ctype + " object should be an integer vector");
        }
        Rcpp::IntegerVector curindices(first);

        // Check for index contents is done inside the fragmented constructor.
        Rcpp::RObject second(inner[1]);
        if (second.sexp_type() != desired_sexp_) {
            auto ctype = get_class_name(seed);
            throw std::runtime_error("second entry of an element of the 'SVT' slot in a " + ctype + " object has an unexpected type");
        }
        InputObject_ curvalues(second);
        size_t nnz = curvalues.size();
        if (nnz != curindices.size()) {
            auto ctype = get_class_name(seed);
            throw std::runtime_error("both vectors of an element of the 'SVT' slot in a " + ctype + " object should have the same length");
        }

        if constexpr(transpose_) {
            auto idx = secondary_extract[c] - 1;
            for (size_t i = 0; i < nnz; ++i) {
                auto ix = curindices[i];
                if (needs_value) {
                    value_ptrs[ix] = curvalues[i];
                }
                if (needs_index) {
                    index_ptrs[ix] = idx;
                }
                ++(counts[ix]);
            }

        } else {
            if (needs_value) {
                auto vptr = value_ptrs[c];
                std::copy(curvalues.begin(), curvalues.end(), vptr);
            }
            if (needs_index) {
                auto iptr = index_ptrs[c];
                for (auto c : curindices) {
                    *iptr = secondary_extract[c];
                    ++iptr;
                }
            }
            counts[c] = nnz;
        }
    }
}

template<bool transpose_, class InputObject_, SEXPTYPE desired_sexp_, typename InputValue_, typename CachedValue_, typename CachedIndex_, typename Index_>
void parse_sparse_matrix(
    Rcpp::RObject seed,
    std::vector<CachedValue_*>& value_ptrs, 
    std::vector<CachedIndex_*>& index_ptrs, 
    std::vector<Index_>& counts, 
    const Rcpp::IntegerVector& secondary_extract)
{
    auto ctype = get_class_name(seed);
    if (ctype != "SVT_SparseMatrix") {
        // Can't be bothered to write a parser for COO_SparseMatrix objects,
        // which are soon-to-be-superceded by SVT_SparseMatrix anyway; so we
        // just forcibly coerce it.
        auto methods_env = Rcpp::Environment::namespace_env("methods");
        Rcpp::Function converter(methods_env["as"]);
        seed = converter(seed, Rcpp::CharacterVector::create("SVT_SparseMatrix"));
    }

    std::string type = Rcpp::as<std::string>(seed.slot("type"));
    if (type == "double") {
        parse_sparse_matrix_internal<transpose_, Rcpp::NumericVector, REALSXP, double>(seed, value_ptrs, index_ptrs, counts, secondary_extract);
    } else if (type == "integer") {
        parse_sparse_matrix_internal<transpose_, Rcpp::IntegerVector, INTSXP, int>(seed, value_ptrs, index_ptrs, counts, secondary_extract);
    } else if (type == "logical") {
        parse_sparse_matrix_internal<transpose_, Rcpp::LogicalVector, LGLSXP, int>(seed, value_ptrs, index_ptrs, counts, secondary_extract);
    } 

    throw std::runtime_error("unsupported type '" + type + "' for a " + ctype);
}

}

#endif
