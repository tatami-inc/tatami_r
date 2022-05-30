#include "Rcpp.h"
#include "raticate/raticate.hpp"

typedef Rcpp::XPtr<raticate::Parsed<double, int> > RatXPtr;

//' @useDynLib raticate.tests
//' @importFrom Rcpp sourceCpp
//[[Rcpp::export(rng=false)]]
SEXP parse(Rcpp::RObject seed) {
    auto out = raticate::convert(seed);
    if (out.matrix == nullptr) {
        throw std::runtime_error("failed to convert R object into a tatami::NumericMatrix");
    }
    return RatXPtr(new raticate::Parsed<double, int>(std::move(out)), true);
}

//[[Rcpp::export(rng=false)]]
int nrow(Rcpp::RObject parsed) {
    RatXPtr ptr(parsed);
    return (ptr->matrix)->nrow();
}

//[[Rcpp::export(rng=false)]]
int ncol(Rcpp::RObject parsed) {
    RatXPtr ptr(parsed);
    return (ptr->matrix)->ncol();
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector row(Rcpp::RObject parsed, int i) {
    RatXPtr ptr(parsed);
    Rcpp::NumericVector output((ptr->matrix)->ncol());
    (ptr->matrix)->row_copy(i, static_cast<double*>(output.begin()));
    return output;
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector column(Rcpp::RObject parsed, int i) {
    RatXPtr ptr(parsed);
    Rcpp::NumericVector output((ptr->matrix)->nrow());
    (ptr->matrix)->column_copy(i, static_cast<double*>(output.begin()));
    return output;
}
