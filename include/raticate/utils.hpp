#ifndef RATICATE_UTILS_HPP
#define RATICATE_UTILS_HPP

/**
 * @file utils.hpp
 *
 * Internal utilities for the **raticate** converters.
 */

#include "Rcpp.h"
#include <string>
#include <utility>
#include <stdexcept>
#include <memory>
#include "tatami/tatami.hpp"

namespace raticate { 

template<typename Data, typename Index>
struct Parsed {
    std::shared_ptr<tatami::Matrix<Data, Index> > matrix;
    Rcpp::List contents;
};

/**
 * @cond
 */
inline std::string make_to_string(const Rcpp::RObject& str) {
    Rcpp::StringVector as_str(str);
    if (as_str.size()!=1) { 
        throw std::runtime_error("input RObject should contain a single string");
    }
    return Rcpp::as<std::string>(as_str[0]);
}

inline Rcpp::RObject get_class_object(const Rcpp::RObject& incoming) {
    if (!incoming.isObject()) {
        throw std::runtime_error("object has no 'class' attribute");
    }
    return incoming.attr("class");
}

inline std::string get_class_name(const Rcpp::RObject& incoming) {
    return make_to_string(get_class_object(incoming));
}

inline std::string extract_class_package(const Rcpp::RObject& classname) {
    if (!classname.hasAttribute("package")) {
        throw std::runtime_error("class name has no 'package' attribute");
    }
    return make_to_string(classname.attr("package"));
}

inline std::pair<std::string, std::string> get_class_package(const Rcpp::RObject& incoming) {
    Rcpp::RObject classname=get_class_object(incoming);
    return std::make_pair(make_to_string(classname), extract_class_package(classname));
}

inline std::string translate_type(int sexp_type) {
    std::string should_be;
    switch(sexp_type) {
        case REALSXP:
            should_be="double";
            break;
        case INTSXP:
            should_be="integer";
            break;
        case LGLSXP:
            should_be="logical";
            break;
        case STRSXP:
            should_be="character";
            break;
        default:
            throw std::runtime_error("unsupported sexptype " + std::to_string(sexp_type));
    }
    return should_be;
}

std::pair<int, int> parse_dims(Rcpp::RObject dims) {
    if (dims.sexp_type()!=INTSXP) {
        throw std::runtime_error("matrix dimensions should be an integer vector");
    }

    Rcpp::IntegerVector d(dims);
    if (d.size()!=2) {
        throw std::runtime_error("matrix dimensions should be of length 2");
    }

    if (d[0]<0 || d[1]<0) {
        throw std::runtime_error("dimensions should be non-negative");
    }

    return std::make_pair(d[0], d[1]);
}
/**
 * @endcond
 */

}

#endif
