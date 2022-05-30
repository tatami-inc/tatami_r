#ifndef RATICATE_HPP
#define RATICATE_HPP

#include "Rcpp.h"
#include "SimpleMatrix.hpp"
#include "SparseArraySeed.hpp"
#include "CSparseMatrix.hpp"
#include "utils.hpp"

namespace raticate {

template<typename Data = double, typename Index = int>
Parsed<Data, Index> convert(Rcpp::RObject x) {
    Parsed<Data, Index> output;

    if (x.isS4()) {
        std::string ctype = get_class_name(x);
        if (ctype == "SparseArraySeed") {
            output = convert_SparseArraySeed<Data, Index>(x);
        } else if (ctype == "dgCMatrix") {
            output = convert_dgCMatrix<Data, Index>(x);
        } else if (ctype == "lgCMatrix") {
            output = convert_lgCMatrix<Data, Index>(x);
        }
    } else if (x.hasAttribute("dim")) {
        output = convert_simple_matrix<Data, Index>(x);
    }

    return output;
}

}

#endif
