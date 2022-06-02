#include "Rcpp.h"
#include "raticate/raticate.hpp"

typedef Rcpp::XPtr<raticate::Parsed<double, int> > RatXPtr;

//' @useDynLib raticate.tests
//' @importFrom Rcpp sourceCpp
//' @export
//[[Rcpp::export(rng=false)]]
SEXP parse(Rcpp::RObject seed) {
    auto out = raticate::parse<double, int>(seed);
    if (out.matrix == nullptr) {
        throw std::runtime_error("failed to parse R object into a tatami::NumericMatrix");
    }
    return RatXPtr(new raticate::Parsed<double, int>(std::move(out)), true);
}

//' @export
//[[Rcpp::export(rng=false)]]
int nrow(Rcpp::RObject parsed) {
    RatXPtr ptr(parsed);
    return (ptr->matrix)->nrow();
}

//' @export
//[[Rcpp::export(rng=false)]]
int ncol(Rcpp::RObject parsed) {
    RatXPtr ptr(parsed);
    return (ptr->matrix)->ncol();
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector row(Rcpp::RObject parsed, int i) {
    RatXPtr ptr(parsed);
    if (i < 1 || i > (ptr->matrix)->nrow()) {
        throw std::runtime_error("requested row index out of range");
    }
    Rcpp::NumericVector output((ptr->matrix)->ncol());
    (ptr->matrix)->row_copy(i - 1, static_cast<double*>(output.begin()));
    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector column(Rcpp::RObject parsed, int i) {
    RatXPtr ptr(parsed);
    if (i < 1 || i > (ptr->matrix)->ncol()) {
        throw std::runtime_error("requested column index out of range");
    }
    Rcpp::NumericVector output((ptr->matrix)->nrow());
    (ptr->matrix)->column_copy(i - 1, static_cast<double*>(output.begin()));
    return output;
}
