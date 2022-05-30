#ifndef RATICATE_HPP
#define RATICATE_HPP

#include "Rcpp.h"
#include "SimpleMatrix.hpp"
#include "SparseArraySeed.hpp"
#include "CSparseMatrix.hpp"
#include "utils.hpp"

namespace raticate {

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
        }
    } else if (x.hasAttribute("dim")) {
        output = parse_simple_matrix<Data, Index>(x);
    }

    return output;
}

}

#endif
