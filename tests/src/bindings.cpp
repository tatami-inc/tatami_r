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
SEXP parse(Rcpp::RObject seed, double cache_size, bool require_min) {
    tatami_r::Options opt;
    opt.maximum_cache_size = cache_size;
    opt.require_minimum_cache = require_min;
    return RatXPtr(new tatami_r::UnknownMatrix<double, int>(seed, opt));
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

//' @export
//[[Rcpp::export(rng=false)]]
bool prefer_row(Rcpp::RObject parsed) {
    RatXPtr ptr(parsed);
    return ptr->prefer_rows();
}

//' @export
//[[Rcpp::export(rng=false)]]
bool is_sparse(Rcpp::RObject parsed) {
    RatXPtr ptr(parsed);
    return ptr->prefer_rows();
}

/******************
 *** Dense full ***
 ******************/

void check_idx(const Rcpp::IntegerVector& idx, int primary) {
    for (auto i : idx) {
        if (i < 1 || i > primary) {
            throw std::runtime_error("requested primary index out of range");
        }
    }
}

std::shared_ptr<tatami::FixedVectorOracle<int> > create_oracle(const Rcpp::IntegerVector& idx) {
    std::vector<int> predictions(idx.begin(), idx.end());
    for (auto& p : predictions) {
        --p;
    }
    return std::make_shared<tatami::FixedVectorOracle<int> >(std::move(predictions));
}

template<bool oracle_, class Extractor_>
Rcpp::NumericVector format_dense_output(Extractor_* ext, int i, int len) {
    Rcpp::NumericVector vec(len);
    auto optr = static_cast<double*>(vec.begin());
    auto iptr = [&]() {
        if constexpr(oracle_) {
            return ext->fetch(optr);
        } else {
            return ext->fetch(i, optr);
        }
    }();
    tatami::copy_n(iptr, len, optr);
    return vec;
}

template<bool sparse_, bool oracle_>
auto create_extractor(const RatXPtr& ptr, bool row, const Rcpp::IntegerVector& idx, bool needs_value = true, bool needs_index = true) {
    tatami::Options opt;
    opt.sparse_extract_value = needs_value;
    opt.sparse_extract_index = needs_index;

    if constexpr(oracle_) {
        return tatami::new_extractor<sparse_, oracle_>(ptr.get(), row, create_oracle(idx), opt);
    } else {
        return tatami::new_extractor<sparse_, oracle_>(ptr.get(), row, false, opt);
    }
}

template<bool oracle_>
Rcpp::List dense_full(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx) {
    RatXPtr ptr(parsed);
    int primary = (row ? ptr->nrow() : ptr->ncol());
    int secondary = (!row ? ptr->nrow() : ptr->ncol());
    check_idx(idx, primary);

    auto ext = create_extractor<false, oracle_>(ptr, row, idx);
    Rcpp::List output(idx.size());
    for (size_t i = 0, end = idx.size(); i < end; ++i) {
        output[i] = format_dense_output<oracle_>(ext.get(), idx[i] - 1, secondary);
    }
    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List myopic_dense_full(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx) {
    return dense_full<false>(std::move(parsed), row, std::move(idx));
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List oracular_dense_full(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx) {
    return dense_full<true>(std::move(parsed), row, std::move(idx));
}

/*******************
 *** Dense block ***
 *******************/

void check_block(int first, int len, int secondary) {
    if (first < 1 || len < 0 || first + len - 1 > secondary) {
        throw std::runtime_error("requested block out of range");
    }
}

template<bool sparse_, bool oracle_>
auto create_extractor(const RatXPtr& ptr, bool row, const Rcpp::IntegerVector& idx, int first, int len, bool needs_value = true, bool needs_index = true) {
    tatami::Options opt;
    opt.sparse_extract_value = needs_value;
    opt.sparse_extract_index = needs_index;

    if constexpr(oracle_) {
        return tatami::new_extractor<sparse_, oracle_>(ptr.get(), row, create_oracle(idx), first - 1, len, opt);
    } else {
        return tatami::new_extractor<sparse_, oracle_>(ptr.get(), row, false, first - 1, len, opt);
    }
}

template<bool oracle_>
Rcpp::List dense_block(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, int first, int len) {
    RatXPtr ptr(parsed);
    int primary = (row ? ptr->nrow() : ptr->ncol());
    int secondary = (!row ? ptr->nrow() : ptr->ncol());
    check_idx(idx, primary);
    check_block(first, len, secondary);

    auto ext = create_extractor<false, oracle_>(ptr, row, idx, first, len);
    Rcpp::List output(idx.size());
    for (size_t i = 0, end = idx.size(); i < end; ++i) {
        output[i] = format_dense_output<oracle_>(ext.get(), idx[i] - 1, len);
    }
    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List myopic_dense_block(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, int first, int len) {
    return dense_block<false>(std::move(parsed), row, std::move(idx), first, len);
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List oracular_dense_block(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, int first, int len) {
    return dense_block<true>(std::move(parsed), row, std::move(idx), first, len);
}

/********************
 *** Dense subset ***
 ********************/

void check_subset(const Rcpp::IntegerVector& subset, int secondary) {
    int last = 0;
    for (auto s : subset) {
        if (s <= last || s > secondary) {
            throw std::runtime_error("subset vector should be sorted and within range");
        }
        last = s;
    }
}

template<bool sparse_, bool oracle_>
auto create_extractor(const RatXPtr& ptr, bool row, const Rcpp::IntegerVector& idx, const Rcpp::IntegerVector& subset, bool needs_value = true, bool needs_index = true) {
    tatami::Options opt;
    opt.sparse_extract_value = needs_value;
    opt.sparse_extract_index = needs_index;

    auto subs = std::make_shared<std::vector<int> >(subset.begin(), subset.end());
    for (auto& s : *subs) {
        --s;
    }

    if constexpr(oracle_) {
        return tatami::new_extractor<sparse_, oracle_>(ptr.get(), row, create_oracle(idx), std::move(subs), opt);
    } else {
        return tatami::new_extractor<sparse_, oracle_>(ptr.get(), row, false, std::move(subs), opt);
    }
}

template<bool oracle_>
Rcpp::List dense_indexed(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, Rcpp::IntegerVector subset) {
    RatXPtr ptr(parsed);
    int primary = (row ? ptr->nrow() : ptr->ncol());
    int secondary = (!row ? ptr->nrow() : ptr->ncol());
    check_idx(idx, primary);
    check_subset(subset, secondary);

    auto ext = create_extractor<false, oracle_>(ptr, row, idx, subset);
    Rcpp::List output(idx.size());
    for (size_t i = 0, end = idx.size(); i < end; ++i) {
        output[i] = format_dense_output<oracle_>(ext.get(), idx[i] - 1, subset.size());
    }

    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List myopic_dense_indexed(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, Rcpp::IntegerVector subset) {
    return dense_indexed<false>(std::move(parsed), row, std::move(idx), std::move(subset));
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List oracular_dense_indexed(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, Rcpp::IntegerVector subset) {
    return dense_indexed<true>(std::move(parsed), row, std::move(idx), std::move(subset));
}

/*******************
 *** Sparse full ***
 *******************/

template<bool oracle_, class Extractor_>
Rcpp::RObject format_sparse_output(Extractor_* ext, int i, double* vbuffer, int* ibuffer, bool needs_value, bool needs_index) {
    if (needs_index && needs_value) {
        auto x = [&](){
            if constexpr(oracle_) {
                return ext->fetch(vbuffer, ibuffer);
            } else {
                return ext->fetch(i, vbuffer, ibuffer);
            }
        }();
        Rcpp::NumericVector outv(x.value, x.value + x.number);
        Rcpp::IntegerVector outi(x.index, x.index + x.number);
        for (auto& i : outi) { ++i; }
        return Rcpp::List::create(Rcpp::Named("index") = outi, Rcpp::Named("value") = outv);

    } else if (needs_index) {
        auto x = [&](){
            if constexpr(oracle_) {
                return ext->fetch(NULL, ibuffer);
            } else {
                return ext->fetch(i, NULL, ibuffer);
            }
        }();
        Rcpp::IntegerVector outi(x.index, x.index + x.number);
        for (auto& i : outi) { ++i; }
        return outi;

    } else if (needs_value) {
        auto x = [&](){
            if constexpr(oracle_) {
                return ext->fetch(vbuffer, NULL);
            } else {
                return ext->fetch(i, vbuffer, NULL);
            }
        }();
        return Rcpp::NumericVector(x.value, x.value + x.number);

    } else {
        auto x = [&](){
            if constexpr(oracle_) {
                return ext->fetch(NULL, NULL);
            } else {
                return ext->fetch(i, NULL, NULL);
            }
        }();
        return Rcpp::IntegerVector::create(x.number);
    }
}

template<bool oracle_>
Rcpp::List sparse_full(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, bool needs_value, bool needs_index) {
    RatXPtr ptr(parsed);
    int primary = (row ? ptr->nrow() : ptr->ncol());
    int secondary = (!row ? ptr->nrow() : ptr->ncol());
    check_idx(idx, primary);

    auto ext = create_extractor<true, oracle_>(ptr, row, idx, needs_value, needs_index);
    Rcpp::List output(idx.size());
    std::vector<double> vbuffer(secondary);
    std::vector<int> ibuffer(secondary);

    for (size_t i = 0, end = idx.size(); i < end; ++i) {
        output[i] = format_sparse_output<oracle_>(ext.get(), idx[i] - 1, vbuffer.data(), ibuffer.data(), needs_value, needs_index);
    }

    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List myopic_sparse_full(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, bool needs_value, bool needs_index) {
    return sparse_full<false>(std::move(parsed), row, std::move(idx), needs_value, needs_index);
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List oracular_sparse_full(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, bool needs_value, bool needs_index) {
    return sparse_full<true>(std::move(parsed), row, std::move(idx), needs_value, needs_index);
}

/********************
 *** Sparse block ***
 ********************/

template<bool oracle_>
Rcpp::List sparse_block(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, int first, int len, bool needs_value, bool needs_index) {
    RatXPtr ptr(parsed);
    int primary = (row ? ptr->nrow() : ptr->ncol());
    int secondary = (!row ? ptr->nrow() : ptr->ncol());
    check_idx(idx, primary);
    check_block(first, len, secondary);

    auto ext = create_extractor<true, oracle_>(ptr, row, idx, first, len, needs_value, needs_index);
    std::vector<double> vbuffer(len);
    std::vector<int> ibuffer(len);
    Rcpp::List output(idx.size());

    for (size_t i = 0, end = idx.size(); i < end; ++i) {
        output[i] = format_sparse_output<oracle_>(ext.get(), idx[i] - 1, vbuffer.data(), ibuffer.data(), needs_value, needs_index);
    }
    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List myopic_sparse_block(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, int first, int len, bool needs_value, bool needs_index) {
    return sparse_block<false>(std::move(parsed), row, std::move(idx), first, len, needs_value, needs_index);
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List oracular_sparse_block(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, int first, int len, bool needs_value, bool needs_index) {
    return sparse_block<true>(std::move(parsed), row, std::move(idx), first, len, needs_value, needs_index);
}

/**********************
 *** Sparse indexed ***
 **********************/

template<bool oracle_>
Rcpp::List sparse_indexed(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, Rcpp::IntegerVector subset, bool needs_value, bool needs_index) {
    RatXPtr ptr(parsed);
    int primary = (row ? ptr->nrow() : ptr->ncol());
    int secondary = (!row ? ptr->nrow() : ptr->ncol());
    check_idx(idx, primary);
    check_subset(subset, secondary);

    auto ext = create_extractor<true, oracle_>(ptr, row, idx, subset, needs_value, needs_index);
    std::vector<double> vbuffer(subset.size());
    std::vector<int> ibuffer(subset.size());
    Rcpp::List output(idx.size());

    for (size_t i = 0, end = idx.size(); i < end; ++i) {
        output[i] = format_sparse_output<oracle_>(ext.get(), idx[i] - 1, vbuffer.data(), ibuffer.data(), needs_value, needs_index);
    }
    return output;
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List myopic_sparse_indexed(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, Rcpp::IntegerVector subset, bool needs_value, bool needs_index) {
    return sparse_indexed<false>(std::move(parsed), row, std::move(idx), std::move(subset), needs_value, needs_index);
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::List oracular_sparse_indexed(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, Rcpp::IntegerVector subset, bool needs_value, bool needs_index) {
    return sparse_indexed<true>(std::move(parsed), row, std::move(idx), std::move(subset), needs_value, needs_index);
}

/****************
 *** Row sums ***
 ****************/

Rcpp::NumericVector collapse_vector(const std::vector<std::vector<double> >& temp) {
    size_t total = 0;
    for (const auto& t : temp) {
        total += t.size();
    }

    Rcpp::NumericVector output(total);
    auto start = output.begin();
    for (const auto& t : temp) {
        std::copy(t.begin(), t.end(), start);
        start += t.size();
    }

    return output;
}

template<bool oracle_>
Rcpp::NumericVector dense_sums(Rcpp::RObject parsed, bool row, int num_threads) {
    RatXPtr ptr(parsed);
    int primary = (row ? ptr->nrow() : ptr->ncol());
    int secondary = (!row ? ptr->nrow() : ptr->ncol());

    std::vector<std::vector<double> > output(num_threads);
    tatami_r::parallelize([&](int w, int start, int len) {
        auto ext = [&]() {
            if constexpr(oracle_) {
                return tatami::new_extractor<false, oracle_>(ptr.get(), row, std::make_shared<tatami::ConsecutiveOracle<int> >(start, len));
            } else {
                return tatami::new_extractor<false, oracle_>(ptr.get(), row, false);
            }
        }();

        std::vector<double> buffer(secondary);
        auto& mine = output[w];
        mine.resize(len);

        for (int i = 0; i < len; ++i) {
            auto iptr = [&]() {
                if constexpr(oracle_) {
                    return ext->fetch(buffer.data());
                } else {
                    return ext->fetch(i + start, buffer.data());
                }
            }();
            mine[i] = std::accumulate(iptr, iptr + secondary, 0.0);
        }
    }, primary, num_threads);

    return collapse_vector(output);
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector myopic_dense_sums(Rcpp::RObject parsed, bool row, int num_threads) {
    return dense_sums<false>(std::move(parsed), row, num_threads);
}

//' @export
//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector oracular_dense_sums(Rcpp::RObject parsed, bool row, int num_threads) {
    return dense_sums<true>(std::move(parsed), row, num_threads);
}

template<bool oracle_>
Rcpp::NumericVector sparse_sums(Rcpp::RObject parsed, bool row, int num_threads) {
    RatXPtr ptr(parsed);
    int primary = (row ? ptr->nrow() : ptr->ncol());
    int secondary = (!row ? ptr->nrow() : ptr->ncol());

    std::vector<std::vector<double> > output(num_threads);
    tatami_r::parallelize([&](int w, int start, int len) {
        auto ext = [&]() {
            if constexpr(oracle_) {
                return tatami::new_extractor<true, oracle_>(ptr.get(), row, std::make_shared<tatami::ConsecutiveOracle<int> >(start, len));
            } else {
                return tatami::new_extractor<true, oracle_>(ptr.get(), row, false);
            }
        }();

        std::vector<double> vbuffer(secondary);
        std::vector<int> ibuffer(secondary);
        auto& mine = output[w];
        mine.resize(len);

        for (int i = 0; i < len; ++i) {
            auto range = [&]() {
                if constexpr(oracle_) {
                    return ext->fetch(vbuffer.data(), ibuffer.data());
                } else {
                    return ext->fetch(i + start, vbuffer.data(), ibuffer.data());
                }
            }();
            mine[i] = std::accumulate(range.value, range.value + range.number, 0.0);
        }
    }, primary, num_threads);

    return collapse_vector(output);
}
