#ifndef RATICATE_HPP
#define RATICATE_HPP

#include "Rcpp.h"
#include "SimpleMatrix.hpp"
#include "SparseArraySeed.hpp"
#include "CSparseMatrix.hpp"
#include "DelayedMatrix.hpp"
#include "DelayedSubset.hpp"
#include "DelayedAperm.hpp"
#include "UnknownMatrix.hpp"
#include "utils.hpp"

/**
 * @file raticate.hpp
 *
 * @brief Parse `Rcpp::RObject`s into **tatami** matrices.
 */

namespace raticate {

/**
 * Parse `Rcpp::RObject`s into **tatami** matrices.
 * Supported matrix types are:
 *
 * - ordinary logical, numeric or integer matrices.
 * - `dgCMatrix` or `lgCMatrix` objects from the **Matrix** package.
 * - `SparseArraySeed` objects from the **DelayedArray** package.
 * - `DelayedMatrix` objects wrapping any of the above, or containing the following delayed operations:
 *    - Subsetting
 *    - Modification of dimnames
 *    - Transposition
 * 
 * @tparam Data Numeric data type for the **tatami** interface.
 * @tparam Index Integer index type for the **tatami** interface.
 * 
 * @param x An R object representing a supported matrix type.
 *
 * @return A `Parsed` object containing a pointer to a parsed `tatami::Matrix`.
 * If parsing was not successful, this pointer will be a `nullptr`.
 */
template<typename Data = double, typename Index = int>
Parsed<Data, Index> parse(Rcpp::RObject x) {
    Parsed<Data, Index> output;

    if (x.isS4()) {
        std::string ctype = get_class_name(x);

        if (ctype == "SparseArraySeed") {
            output = parse_SparseArraySeed<Data, Index>(x);
        } else if (ctype == "dgCMatrix") {
            output = parse_dgCMatrix<Data, Index>(x);
        } else if (ctype == "lgCMatrix") {
            output = parse_lgCMatrix<Data, Index>(x);

        } else if (ctype == "DelayedMatrix") {
            output = parse_DelayedMatrix<Data, Index>(x);
        } else if (ctype == "DelayedSetDimnames") {
            output = parse_DelayedMatrix<Data, Index>(x); // just forward onto the seed.
        } else if (ctype == "DelayedSubset") {
            output = parse_DelayedSubset<Data, Index>(x);
        } else if (ctype == "DelayedAperm") {
            output = parse_DelayedAperm<Data, Index>(x);
        }

    } else if (x.hasAttribute("dim")) {
        output = parse_simple_matrix<Data, Index>(x);
    }

    if (output.matrix == nullptr) {
        // No need to set contents here, as the matrix itself holds the Rcpp::RObject.
        output.matrix.reset(new UnknownMatrix<Data, Index>(x));
    }

    return output;
}

}

#endif
