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
    std::vector<Index_>& counts)
{
    Rcpp::RObject raw_svt = seed.slot("SVT");
    if (raw_svt == R_NilValue) {
        return;
    }

    Rcpp::List svt(raw_svt);
    int NC = svt.size();
    bool needs_value = !value_ptrs.empty();
    bool needs_index = !index_ptrs.empty();

    // Note that non-empty value_ptrs and index_ptrs may be longer than the
    // number of rows/columns in the SVT matrix, due to the reuse of slabs.

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
        if (nnz != static_cast<size_t>(curindices.size())) {
            auto ctype = get_class_name(seed);
            throw std::runtime_error("both vectors of an element of the 'SVT' slot in a " + ctype + " object should have the same length");
        }

        if constexpr(transpose_) {
            for (size_t i = 0; i < nnz; ++i) {
                auto ix = curindices[i];
                auto& shift = counts[ix];
                if (needs_value) {
                    value_ptrs[ix][shift] = curvalues[i];
                }
                if (needs_index) {
                    index_ptrs[ix][shift] = c;
                }
                ++shift;
            }

        } else {
            if (needs_value) {
                std::copy(curvalues.begin(), curvalues.end(), value_ptrs[c]);
            }
            if (needs_index) {
                std::copy(curindices.begin(), curindices.end(), index_ptrs[c]);
            }
            counts[c] = nnz;
        }
    }
}

template<bool transpose_, typename CachedValue_, typename CachedIndex_, typename Index_>
void parse_sparse_matrix(
    Rcpp::RObject seed,
    std::vector<CachedValue_*>& value_ptrs, 
    std::vector<CachedIndex_*>& index_ptrs, 
    std::vector<Index_>& counts)
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
        parse_sparse_matrix_internal<transpose_, Rcpp::NumericVector, REALSXP, double>(seed, value_ptrs, index_ptrs, counts);
    } else if (type == "integer") {
        parse_sparse_matrix_internal<transpose_, Rcpp::IntegerVector, INTSXP, int>(seed, value_ptrs, index_ptrs, counts);
    } else if (type == "logical") {
        parse_sparse_matrix_internal<transpose_, Rcpp::LogicalVector, LGLSXP, int>(seed, value_ptrs, index_ptrs, counts);
    } else {
        throw std::runtime_error("unsupported type '" + type + "' for a " + ctype);
    }
}

}

#endif
