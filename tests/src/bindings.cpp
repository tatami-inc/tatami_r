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
Rcpp::NumericVector row_subset(Rcpp::RObject parsed, int i, int first, int last) {
    RatXPtr ptr(parsed);
    if (i < 1 || i > (ptr->matrix)->nrow()) {
        throw std::runtime_error("requested row index out of range");
    }
    if (first < 1 || first > last || last > (ptr->matrix)->ncol()) {
        throw std::runtime_error("requested subset indices out of range");
    }
    Rcpp::NumericVector output(last - first + 1);
    (ptr->matrix)->row_copy(i - 1, static_cast<double*>(output.begin()), first - 1, last);
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

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector column_subset(Rcpp::RObject parsed, int i, int first, int last) {
    RatXPtr ptr(parsed);
    if (i < 1 || i > (ptr->matrix)->ncol()) {
        throw std::runtime_error("requested column index out of range");
    }
    if (first < 1 || first > last || last > (ptr->matrix)->nrow()) {
        throw std::runtime_error("requested subset indices out of range");
    }
    Rcpp::NumericVector output(last - first + 1);
    (ptr->matrix)->column_copy(i - 1, static_cast<double*>(output.begin()), first - 1, last);
    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List rows(Rcpp::RObject parsed) {
    RatXPtr ptr(parsed);

    auto wrk = (ptr->matrix)->new_workspace(true);
    size_t nr = (ptr->matrix)->nrow();
    size_t nc = (ptr->matrix)->ncol();
    Rcpp::List output(nr);

    for (size_t r = 0; r < nr; ++r) {
        Rcpp::NumericVector current(nc);
        (ptr->matrix)->row_copy(r, static_cast<double*>(current.begin()), wrk.get());
        output[r] = current;
    }

    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List rows_subset(Rcpp::RObject parsed, int first, int last) {
    RatXPtr ptr(parsed);

    auto wrk = (ptr->matrix)->new_workspace(true);
    size_t nr = (ptr->matrix)->nrow();
    size_t len = last - first + 1;
    Rcpp::List output(nr);

    for (size_t r = 0; r < nr; ++r) {
        Rcpp::NumericVector current(len);
        (ptr->matrix)->row_copy(r, static_cast<double*>(current.begin()), first - 1, last, wrk.get());
        output[r] = current;
    }

    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List columns(Rcpp::RObject parsed) {
    RatXPtr ptr(parsed);

    auto wrk = (ptr->matrix)->new_workspace(false);
    size_t nr = (ptr->matrix)->nrow();
    size_t nc = (ptr->matrix)->ncol();
    Rcpp::List output(nc);

    for (size_t c = 0; c < nc; ++c) {
        Rcpp::NumericVector current(nr);
        (ptr->matrix)->column_copy(c, static_cast<double*>(current.begin()), wrk.get());
        output[c] = current;
    }

    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List columns_subset(Rcpp::RObject parsed, int first, int last) {
    RatXPtr ptr(parsed);

    auto wrk = (ptr->matrix)->new_workspace(false);
    size_t nc = (ptr->matrix)->ncol();
    size_t len = last - first + 1;
    Rcpp::List output(nc);

    for (size_t c = 0; c < nc; ++c) {
        Rcpp::NumericVector current(len);
        (ptr->matrix)->column_copy(c, static_cast<double*>(current.begin()), first - 1, last, wrk.get());
        output[c] = current;
    }

    return output;
}


