#include "Rcpp.h"
#include <vector>
#include <algorithm>
#include <iostream>

#ifdef TEST_CUSTOM_PARALLEL
#define TATAMI_R_PARALLELIZE_UNKNOWN
#include "tatami_r/parallelize.hpp"
#define TATAMI_CUSTOM_PARALLEL tatami_r::parallelize
#endif

#include "tatami_r/tatami_r.hpp"

typedef Rcpp::XPtr<tatami::Matrix<double, int> > RatXPtr;

//' @useDynLib raticate.tests
//' @importFrom Rcpp sourceCpp
//' @export
//[[Rcpp::export(rng=false)]]
SEXP parse(Rcpp::RObject seed) {
    return RatXPtr(new tatami_r::UnknownMatrix<double, int>(seed));
}

//' @export
//[[Rcpp::export(rng=false)]]
int nrow(Rcpp::RObject parsed) {
    RatXPtr ptr(parsed);
    return ptr->nrow();
}

//' @export
//[[Rcpp::export(rng=false)]]
int ncol(Rcpp::RObject parsed) {
    RatXPtr ptr(parsed);
    return ptr->ncol();
}

/******************************************
 *** Dense single row/column extractors ***
 ******************************************/

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector row(Rcpp::RObject parsed, int i) {
    RatXPtr ptr(parsed);
    if (i < 1 || i > ptr->nrow()) {
        throw std::runtime_error("requested row index out of range");
    }
    Rcpp::NumericVector output(ptr->ncol());
    ptr->dense_row()->fetch(i - 1, static_cast<double*>(output.begin()));
    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector row_subset(Rcpp::RObject parsed, int i, int first, int last) {
    RatXPtr ptr(parsed);
    if (i < 1 || i > ptr->nrow()) {
        throw std::runtime_error("requested row index out of range");
    }
    if (first < 1 || first > last || last > ptr->ncol()) {
        throw std::runtime_error("requested subset indices out of range");
    }
    Rcpp::NumericVector output(last - first + 1);
    ptr->dense_row(first - 1, output.size())->fetch_copy(i - 1, static_cast<double*>(output.begin()));
    return output;
}


//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector column(Rcpp::RObject parsed, int i) {
    RatXPtr ptr(parsed);
    if (i < 1 || i > ptr->ncol()) {
        throw std::runtime_error("requested column index out of range");
    }
    Rcpp::NumericVector output(ptr->nrow());
    ptr->dense_column()->fetch_copy(i - 1, static_cast<double*>(output.begin()));
    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector column_subset(Rcpp::RObject parsed, int i, int first, int last) {
    RatXPtr ptr(parsed);
    if (i < 1 || i > ptr->ncol()) {
        throw std::runtime_error("requested column index out of range");
    }
    if (first < 1 || first > last || last > ptr->nrow()) {
        throw std::runtime_error("requested subset indices out of range");
    }
    Rcpp::NumericVector output(last - first + 1);
    ptr->dense_column(first - 1, output.size())->fetch_copy(i - 1, static_cast<double*>(output.begin()));
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
    if (i < 1 || i > ptr->nrow()) {
        throw std::runtime_error("requested row index out of range");
    }
    auto output = ptr->sparse_row()->fetch(i - 1);
    return format_sparse_range(output);
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List sparse_row_subset(Rcpp::RObject parsed, int i, int first, int last) {
    RatXPtr ptr(parsed);
    if (i < 1 || i > ptr->nrow()) {
        throw std::runtime_error("requested row index out of range");
    }
    if (first < 1 || first > last || last > ptr->ncol()) {
        throw std::runtime_error("requested subset indices out of range");
    }
    auto output = ptr->sparse_row(first - 1, last - first + 1)->fetch(i - 1);
    return format_sparse_range(output);
}


//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List sparse_column(Rcpp::RObject parsed, int i) {
    RatXPtr ptr(parsed);
    if (i < 1 || i > ptr->ncol()) {
        throw std::runtime_error("requested column index out of range");
    }
    auto output = ptr->sparse_column()->fetch(i - 1);
    return format_sparse_range(output);
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List sparse_column_subset(Rcpp::RObject parsed, int i, int first, int last) {
    RatXPtr ptr(parsed);
    if (i < 1 || i > ptr->ncol()) {
        throw std::runtime_error("requested column index out of range");
    }
    if (first < 1 || first > last || last > ptr->nrow()) {
        throw std::runtime_error("requested subset indices out of range");
    }
    auto output = ptr->sparse_column(first - 1, last - first + 1)->fetch(i - 1);
    return format_sparse_range(output);
}

/********************************************
 *** Dense multiple row/column extractors ***
 ********************************************/

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List rows(Rcpp::RObject parsed) {
    RatXPtr ptr(parsed);

    size_t nr = ptr->nrow();
    size_t nc = ptr->ncol();
    Rcpp::List output(nr);

    auto wrk = ptr->dense_row();
    for (size_t r = 0; r < nr; ++r) {
        Rcpp::NumericVector current(nc);
        wrk->fetch_copy(r, static_cast<double*>(current.begin()));
        output[r] = current;
    }

    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List rows_subset(Rcpp::RObject parsed, int first, int last) {
    RatXPtr ptr(parsed);

    size_t nr = ptr->nrow();
    size_t len = last - first + 1;
    Rcpp::List output(nr);

    auto wrk = ptr->dense_row(first - 1, last - first + 1);
    for (size_t r = 0; r < nr; ++r) {
        Rcpp::NumericVector current(len);
        wrk->fetch_copy(r, static_cast<double*>(current.begin()));
        output[r] = current;
    }

    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List columns(Rcpp::RObject parsed) {
    RatXPtr ptr(parsed);

    size_t nr = ptr->nrow();
    size_t nc = ptr->ncol();
    Rcpp::List output(nc);

    auto wrk = ptr->dense_column();
    for (size_t c = 0; c < nc; ++c) {
        Rcpp::NumericVector current(nr);
        wrk->fetch_copy(c, static_cast<double*>(current.begin()));
        output[c] = current;
    }

    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List columns_subset(Rcpp::RObject parsed, int first, int last) {
    RatXPtr ptr(parsed);

    size_t nc = ptr->ncol();
    size_t len = last - first + 1;
    Rcpp::List output(nc);

    auto wrk = ptr->dense_column(first - 1, last - first + 1);
    for (size_t c = 0; c < nc; ++c) {
        Rcpp::NumericVector current(len);
        wrk->fetch_copy(c, static_cast<double*>(current.begin()));
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

    size_t nr = ptr->nrow();
    Rcpp::List output(nr);

    auto wrk = ptr->sparse_row();
    for (size_t r = 0; r < nr; ++r) {
        auto current = wrk->fetch(r);
        output[r] = format_sparse_range(current);
    }

    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List sparse_rows_subset(Rcpp::RObject parsed, int first, int last) {
    RatXPtr ptr(parsed);

    size_t nr = ptr->nrow();
    Rcpp::List output(nr);

    auto wrk = ptr->sparse_row(first - 1, last - first + 1);
    for (size_t r = 0; r < nr; ++r) {
        auto current = wrk->fetch(r);
        output[r] = format_sparse_range(current);
    }

    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List sparse_columns(Rcpp::RObject parsed) {
    RatXPtr ptr(parsed);

    size_t nc = ptr->ncol();
    Rcpp::List output(nc);

    auto wrk = ptr->sparse_column();
    for (size_t c = 0; c < nc; ++c) {
        auto current = wrk->fetch(c);
        output[c] = format_sparse_range(current);
    }

    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List sparse_columns_subset(Rcpp::RObject parsed, int first, int last) {
    RatXPtr ptr(parsed);

    size_t nc = ptr->ncol();
    Rcpp::List output(nc);

    auto wrk = ptr->sparse_column(first - 1, last - first + 1);
    for (size_t c = 0; c < nc; ++c) {
        auto current = wrk->fetch(c);
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
    auto out = tatami::row_sums(ptr.get());
    return Rcpp::NumericVector(out.begin(), out.end());
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector rowsums_manual(Rcpp::RObject parsed) {
    RatXPtr ptr(parsed);
    size_t NR = ptr->nrow();
    std::vector<double> output(NR);

#ifdef TEST_CUSTOM_PARALLEL
    int nthreads = 3;
#else
    int nthreads = 1;
#endif

    tatami::parallelize([&](int, int start, int length) -> void {
        auto wrk = ptr->dense_row();
        for (size_t r = start, e = start + length; r < e; ++r) {
            auto current = wrk->fetch(r);
            output[r] = std::accumulate(current.begin(), current.end(), 0.0);
        }
    }, NR, nthreads); 

    return Rcpp::NumericVector(output.begin(), output.end());
}

/*************************
 *** Guided extractors ***
 *************************/

template<class Extractor_>
std::vector<int> prepare_indices(const Rcpp::IntegerVector& targets, Extractor_& wrk) {
    std::vector<int> copy(targets.begin(), targets.end());
    for (auto& x : copy) { --x; }
    wrk->set_oracle(std::make_unique<tatami::FixedOracle<int> >(copy.data(), copy.size()));
    return copy;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List dense_rows_guided(Rcpp::RObject parsed, Rcpp::IntegerVector targets) {
    RatXPtr ptr(parsed);
    auto wrk = ptr->dense_row();
    auto copy = prepare_indices(targets, wrk);

    size_t nc = ptr->ncol();
    Rcpp::List output(copy.size());
    size_t counter = 0;
    for (auto r : copy) {
        Rcpp::NumericVector current(nc);
        wrk->fetch_copy(r, static_cast<double*>(current.begin()));
        output[counter] = current;
        ++counter;
    }

    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List dense_columns_guided(Rcpp::RObject parsed, Rcpp::IntegerVector targets) {
    RatXPtr ptr(parsed);
    auto wrk = ptr->dense_column();
    auto copy = prepare_indices(targets, wrk);

    size_t nr = ptr->nrow();
    size_t counter = 0;
    Rcpp::List output(copy.size());
    for (auto c : copy) {
        Rcpp::NumericVector current(nr);
        wrk->fetch_copy(c, static_cast<double*>(current.begin()));
        output[counter] = current;
        ++counter;
    }

    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List sparse_rows_guided(Rcpp::RObject parsed, Rcpp::IntegerVector targets) {
    RatXPtr ptr(parsed);
    auto wrk = ptr->sparse_row();
    auto copy = prepare_indices(targets, wrk);

    Rcpp::List output(copy.size());
    size_t counter = 0;
    for (auto r : copy) {
        output[counter] = format_sparse_range(wrk->fetch(r));
        ++counter;
    }

    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List sparse_columns_guided(Rcpp::RObject parsed, Rcpp::IntegerVector targets) {
    RatXPtr ptr(parsed);
    auto wrk = ptr->sparse_column();
    auto copy = prepare_indices(targets, wrk);

    size_t counter = 0;
    Rcpp::List output(copy.size());
    for (auto c : copy) {
        output[counter] = format_sparse_range(wrk->fetch(c));
        ++counter;
    }

    return output;
}
