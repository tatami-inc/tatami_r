#include "Rcpp.h"
#include <vector>
#include <algorithm>
#include <iostream>

#ifdef TEST_CUSTOM_PARALLEL
#define RATICATE_PARALLELIZE_UNKNOWN

template<class Function> 
void run(size_t n, Function f);
#define TATAMI_CUSTOM_PARALLEL run
#endif

#include "raticate/raticate.hpp"

#ifdef TEST_CUSTOM_PARALLEL
template<class Function> 
void run(size_t n, Function f) {
    raticate::parallel_coordinator().run<double, int>(n, f, 3);
}
#endif

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

/******************************************
 *** Dense single row/column extractors ***
 ******************************************/

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

/*******************************************
 *** Sparse single row/column extractors ***
 *******************************************/

template<class RangeCopy>
Rcpp::List format_sparse_range(const RangeCopy& x) {
    Rcpp::IntegerVector idx(x.index.begin(), x.index.end());
    for (auto& i : idx) { ++i; }
    return Rcpp::List::create(
        Rcpp::Named("index") = idx,
        Rcpp::Named("value") = Rcpp::NumericVector(x.value.begin(), x.value.end())
    );
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List sparse_row(Rcpp::RObject parsed, int i) {
    RatXPtr ptr(parsed);
    if (i < 1 || i > (ptr->matrix)->nrow()) {
        throw std::runtime_error("requested row index out of range");
    }
    auto output = (ptr->matrix)->sparse_row(i - 1);
    return format_sparse_range(output);
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List sparse_row_subset(Rcpp::RObject parsed, int i, int first, int last) {
    RatXPtr ptr(parsed);
    if (i < 1 || i > (ptr->matrix)->nrow()) {
        throw std::runtime_error("requested row index out of range");
    }
    if (first < 1 || first > last || last > (ptr->matrix)->ncol()) {
        throw std::runtime_error("requested subset indices out of range");
    }
    auto output = (ptr->matrix)->sparse_row(i - 1, first - 1, last);
    return format_sparse_range(output);
}


//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List sparse_column(Rcpp::RObject parsed, int i) {
    RatXPtr ptr(parsed);
    if (i < 1 || i > (ptr->matrix)->ncol()) {
        throw std::runtime_error("requested column index out of range");
    }
    auto output = (ptr->matrix)->sparse_column(i - 1);
    return format_sparse_range(output);
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List sparse_column_subset(Rcpp::RObject parsed, int i, int first, int last) {
    RatXPtr ptr(parsed);
    if (i < 1 || i > (ptr->matrix)->ncol()) {
        throw std::runtime_error("requested column index out of range");
    }
    if (first < 1 || first > last || last > (ptr->matrix)->nrow()) {
        throw std::runtime_error("requested subset indices out of range");
    }
    auto output = (ptr->matrix)->sparse_column(i - 1, first - 1, last);
    return format_sparse_range(output);
}

/********************************************
 *** Dense multiple row/column extractors ***
 ********************************************/

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

/********************************************
 *** Dense multiple row/column extractors ***
 ********************************************/

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List sparse_rows(Rcpp::RObject parsed) {
    RatXPtr ptr(parsed);

    auto wrk = (ptr->matrix)->new_workspace(true);
    size_t nr = (ptr->matrix)->nrow();
    size_t nc = (ptr->matrix)->ncol();
    Rcpp::List output(nr);

    for (size_t r = 0; r < nr; ++r) {
        auto current = (ptr->matrix)->sparse_row(r, wrk.get());
        output[r] = format_sparse_range(current);
    }

    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List sparse_rows_subset(Rcpp::RObject parsed, int first, int last) {
    RatXPtr ptr(parsed);

    auto wrk = (ptr->matrix)->new_workspace(true);
    size_t nr = (ptr->matrix)->nrow();
    size_t len = last - first + 1;
    Rcpp::List output(nr);

    for (size_t r = 0; r < nr; ++r) {
        auto current = (ptr->matrix)->sparse_row(r, first - 1, last, wrk.get());
        output[r] = format_sparse_range(current);
    }

    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List sparse_columns(Rcpp::RObject parsed) {
    RatXPtr ptr(parsed);

    auto wrk = (ptr->matrix)->new_workspace(false);
    size_t nr = (ptr->matrix)->nrow();
    size_t nc = (ptr->matrix)->ncol();
    Rcpp::List output(nc);

    for (size_t c = 0; c < nc; ++c) {
        auto current = (ptr->matrix)->sparse_column(c, wrk.get());
        output[c] = format_sparse_range(current);
    }

    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List sparse_columns_subset(Rcpp::RObject parsed, int first, int last) {
    RatXPtr ptr(parsed);

    auto wrk = (ptr->matrix)->new_workspace(false);
    size_t nc = (ptr->matrix)->ncol();
    size_t len = last - first + 1;
    Rcpp::List output(nc);

    for (size_t c = 0; c < nc; ++c) {
        auto current = (ptr->matrix)->sparse_column(c, first - 1, last, wrk.get());
        output[c] = format_sparse_range(current);
    }

    return output;
}

/*********************************
 *** Parallelizable extractors ***
 *********************************/

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector rowsums(Rcpp::RObject parsed) {
    RatXPtr ptr(parsed);
    auto out = tatami::row_sums((ptr->matrix).get());
    return Rcpp::NumericVector(out.begin(), out.end());
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector rowsums_manual(Rcpp::RObject parsed) {
    RatXPtr ptr(parsed);
    size_t NR = (ptr->matrix)->nrow();
    std::vector<double> output(NR);

#ifndef TATAMI_CUSTOM_PARALLEL
    for (size_t r = 0; r < NR; ++r) {
#else
    TATAMI_CUSTOM_PARALLEL(NR, [&](size_t first, size_t last) -> void {
    for (size_t r = first; r < last; ++r) {
#endif

        auto current = (ptr->matrix)->row(r);
        output[r] = std::accumulate(current.begin(), current.end(), 0.0);

#ifndef TATAMI_CUSTOM_PARALLEL
    }
#else
    }});
#endif

    return Rcpp::NumericVector(output.begin(), output.end());
}
