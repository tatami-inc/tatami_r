#ifndef RATICATE_HPP
#define RATICATE_HPP

#include "Rcpp.h"
#include "SimpleMatrix.hpp"
#include "SparseArraySeed.hpp"
#include "CSparseMatrix.hpp"
#include "DelayedMatrix.hpp"
#include "DelayedSubset.hpp"
#include "DelayedAperm.hpp"
#include "DelayedAbind.hpp"
#include "UnknownMatrix.hpp"
#include "utils.hpp"
#include "parallelize.hpp"

/**
 * @file raticate.hpp
 *
 * @brief Parse `Rcpp::RObject`s into **tatami** matrices.
 */

namespace raticate {

/**
 * Parse `Rcpp::RObject`s into **tatami** matrices.
 * Natively supported matrix types are:
 *
 * - ordinary logical, numeric or integer matrices.
 * - `dgCMatrix` or `lgCMatrix` objects from the **Matrix** package.
 * - `SparseArraySeed` objects from the [**DelayedArray**](https://github.com/Bioconductor/DelayedArray) package.
 * - `DelayedMatrix` objects wrapping any of the above, or containing the following delayed operations:
 *    - Subsetting (as a `DelayedSubset` instance)
 *    - Modification of dimnames (as a `DelayedSetDimnames` instance)
 *    - Transposition (as a `DelayedAperm` instance)
 *    - Combining (as a `DelayedAbind` instance)
 * 
 * For all other objects, we call `DelayedArray::extract_array()` to extract an appropriate slice of the matrix.
 * This is quite a bit slower as it involves a call into the R runtime.
 * 
 * @tparam Data Numeric data type for the **tatami** interface, typically `double`.
 * @tparam Index Integer index type for the **tatami** interface, typically `int`.
 * 
 * @param x An R object representing a supported matrix type.
 *
 * @return A `Parsed` object containing a pointer to a parsed `tatami::Matrix`.
 * If parsing was not successful, this pointer will be a `nullptr`.
 */
template<typename Data, typename Index>
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
        } else if (ctype == "DelayedAbind") {
            output = parse_DelayedAbind<Data, Index>(x);
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