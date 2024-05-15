#ifndef TATAMI_R_DENSE_EXTRACTOR_HPP
#define TATAMI_R_DENSE_EXTRACTOR_HPP

#include "Rcpp.h"
#include "tatami/tatami.hpp"
#include "tatami_chunked/tatami_chunked.hpp"
#include "dense_matrix.hpp"

#include <vector>
#include <stdexcept>
#include <type_traits>

namespace tatami_r {

namespace UnknownMatrix_internal {

/********************
 *** Core classes ***
 ********************/

template<bool oracle_, typename Index_> 
struct SoloDenseCore {
    SoloDenseBase(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& dense_extractor,
        tatami::MaybeOracle<oracle_, Index_> ora,
        Rcpp::IntegerVector non_target_extract, 
        bool by_column,
        [[maybe_unused]] const std::vector<Index_>& ticks,
        [[maybe_unused]] const std::vector<Index_>& map,
        [[maybe_unused]] const tatami_chunked::SlabCacheStats& stats) :
        mat(mat),
        dense_extractor(dense_extractor),
        extract_args(2),
        by_column(by_column)
    {
        extract_args[static_cast<int>(!by_column)] = non_target_extract;
    }

private:
    const Rcpp::RObject& mat;
    const Rcpp::Function& dense_extractor;
    Rcpp::List extract_args;
    bool by_column;

public:
    template<typename Value_>
    void fetch_raw(Index_ i, Value_* buffer) {
        if constexpr(oracle_) {
            i = cache.cache.next();
        }

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
        // This involves some Rcpp initializations, so we lock it just in case.
        auto& mexec = executor();
        mexec.run([&]() -> void {
#endif

        extract_args[static_cast<int>(by_column)] = Rcpp::IntegerVector::create(i + 1);
        auto obj = dense_extractor(mat, extract_args);
        if (by_column) {
            parse_dense_matrix<false>(obj, buffer, 0, 0, non_target_length, 1);
        } else {
            parse_dense_matrix<true>(obj, buffer, 0, 0, 1, non_target_length);
        }

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
        });
#endif
    }
};

template<typename Index_, typename CachedValue_>
struct MyopicDenseCore {
    MyopicDenseBase(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& dense_extractor,
        [[maybe_unused]] tatami::MaybeOracle<oracle_, Index_> ora,
        Rcpp::IntegerVector non_target_extract, 
        bool by_column,
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats& stats) :
        mat(mat),
        dense_extractor(dense_extractor),
        extract_args(2),
        by_column(by_column),
        chunk_ticks(ticks),
        chunk_map(map),
        non_target_length(non_target_extract.size()),
        factory(stats),
        cache(max_slabs)
    {
        extract_args[static_cast<int>(!by_column)] = non_target_extract;
    }

private:
    const Rcpp::RObject& mat;
    const Rcpp::Function& dense_extractor;
    Rcpp::List extract_args;
    bool by_column;

    const std::vector<Index_>& chunk_ticks;
    const std::vector<Index_>& chunk_map;
    size_t non_target_length;

    tatami_chunked::DenseSlabFactory<CachedValue_> factory;
    typedef typename decltype(factory)::Slab;
    tatami_chunked::LruSlabCache<Index_, Slab> cache;

public:
    template<typename Value_>
    void fetch_raw(Index_ i, Value_* buffer) {
        auto chosen = chunk_map[i];

        const auto& slab = cache.find(
            chosen,
            [&]() -> Slab {
                return factory.create();
            },
            [&](Index_ id, Slab& cache) {
#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
                // This involves some Rcpp initializations, so we lock it just in case.
                auto& mexec = executor();
                mexec.run([&]() -> void {
#endif

                auto chunk_start = chunk_ticks[id];
                size_t chunk_len = chunk_ticks[id + 1] - chunk_start;
                Rcpp::IntegerVector primary_extract(chunk_len);
                std::iota(primary_extract.begin(), primary_extract.end(), chunk_start + 1);
                extract_args[static_cast<int>(by_column)] = primary_extract;
                auto obj = dense_extractor(mat, extract_args);
                if (by_column) {
                    parse_dense_matrix<false>(obj, cache.data, 0, 0, non_target_length, chunk_len);
                } else {
                    parse_dense_matrix<true>(obj, cache.data, 0, 0, chunk_len, non_target_length);
                }

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
                });
#endif
            }
        );

        auto src = slab.data + static_cast<size_t>(i - chunk_ticks[chosen]) * non_target_length;
        std::copy_n(src, non_target_length, buffer);
    }
};

template<typename Index_, typename CachedValue_>
struct OracularDenseCore {
    OracularDenseBase(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& dense_extractor,
        tatami::MaybeOracle<oracle_, Index_> ora,
        Rcpp::IntegerVector non_target_extract, 
        bool by_column,
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats& stats) :
        mat(mat),
        dense_extractor(dense_extractor),
        extract_args(2),
        by_column(by_column),
        chunk_ticks(ticks),
        chunk_map(map),
        non_target_length(non_target_extract.size()),
        factory(stats),
        cache(std::move(ora), max_slabs)
    {
        extract_args[static_cast<int>(!by_column)] = non_target_extract;
    }

private:
    const Rcpp::RObject& mat;
    const Rcpp::Function& dense_extractor;
    Rcpp::List extract_args;
    bool by_column;

    const std::vector<Index_>& chunk_ticks;
    const std::vector<Index_>& chunk_map;
    size_t non_target_length;

    tatami_chunked::DenseSlabFactory<CachedValue_> factory;
    typedef typename decltype(factory)::Slab;
    tatami_chunked::OracularSlabCache<Index_, Index_, Slab> cache;

public
    template<typename Value_>
    void fetch_raw(Index_, Value_* buffer) {
        auto res = cache.cache.next(
            [&](Index_ i) -> std::pair<Index_, Index_> {
                auto chosen = chunk_map[i];
                return std::make_pair(chosen, static_cast<Index_>(i - chunk_ticks[chosen]));
            },
            [&]() -> Slab {
                return factory.create();
            },
            [&](std::vector<std::pair<Index_, Slab*> >& to_populate) {
                // Sorting them so that the indices are in order.
                if (!std::is_sorted(to_populate.begin(), to_populate.end(), [&](const std::pair<Index_, Slab*>& left, const std::pair<Index_, Slab*> right) { return left.first < right.first; })) {
                    std::sort(to_populate.begin(), to_populate.end(), [&](const std::pair<Index_, Slab*>& left, const std::pair<Index_, Slab*> right) { return left.first < right.first; });
                }

                Index_ total_len = 0;
                for (const auto& p : to_populate) {
                    total_len += chunk_ticks[p.first + 1] - chunk_ticks[p.first];
                }

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
                // This involves some Rcpp initializations, so we lock it just in case.
                auto& mexec = executor();
                mexec.run([&]() -> void {
#endif

                Rcpp::IntegerVector primary_extract(total_len);
                Index_ current = 0;
                for (const auto& p : to_populate) {
                    Index_ chunk_start = chunk_ticks[p.first];
                    Index_ chunk_len = chunk_ticks[p.first + 1] - chunk_start;
                    auto start = primary_extract.begin() + current;
                    std::iota(start, start + chunk_len, chunk_ticks[p.first] + 1);
                    current += chunk_len;
                }

                extract_args[static_cast<int>(by_column)] = primary_extract;
                auto obj = dense_extractor(mat, extract_args);

                current = 0;
                for (const auto& p : to_populate) {
                    auto chunk_start = chunk_ticks[p.first];
                    Index_ chunk_len = chunk_ticks[p.first + 1] - chunk_start;
                    if (by_column) {
                        parse_dense_matrix<false>(obj, p.second->data, 0, current, non_target_length, chunk_len);
                    } else {
                        parse_dense_matrix<true>(obj, p.second->data, current, 0, chunk_len, non_target_length);
                    }
                    current += chunk_len;
                }

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
                });
#endif
            }
        );

        size_t shift = non_target_length * static_cast<size_t>(res.second); // cast to size_t to avoid overflow.
        std::copy_n(res.first->data + shift, non_target_length, buffer);
    }
};

template<bool solo_ bool oracle_, typename Index_, typename CachedValue_>
using DenseCore = typename std::conditional<solo_,
    SoloDenseCore<oracle_, Index_>,
    typename std::conditional<oracle_,
        OracularDenseCore<Index_, CachedValue_>,
        MyopicDenseCore<Index_, CachedValue_>
    >::type
>::type;

/*************************
 *** Extractor classes ***
 *************************/

template<bool solo_, bool oracle_, typename Value_, typename Index_, typename CachedValue_>
struct DenseFull : public tatami::DenseExtractor<oracle_, Value_, Index_> {
    DenseFull(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& dense_extractor,
        tatami::MaybeOracle<oracle_, Index_> ora,
        Index_ non_target_dim,
        bool by_column,
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats& stats) :
        core(
            mat,
            dense_extractor,
            std::move(ora),
            [&]() {
                Rcpp::IntegerVector output(non_target_dim);
                std::iota(output.begin(), output.end(), 1);
                return output;
            }(),
            by_column,
            ticks,
            map,
            stats
        )
    {}

private:
    DenseCore<solo_, oracle_, Index_, CachedValue_> core;

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        core->fetch_raw(i, buffer);
        return buffer;
    }
};

template<bool solo_, bool oracle_, typename Value_, typename Index_, typename CachedValue_>
struct DenseBlock : public tatami::DenseExtractor<oracle_, Value_, Index_> {
    DenseBlock(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& dense_extractor,
        tatami::MaybeOracle<oracle_, Index_> ora,
        Index_ block_start,
        Index_ block_length,
        bool by_column,
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats& stats) :
        core(
            mat,
            dense_extractor,
            std::move(ora),
            [&]() {
                Rcpp::IntegerVector output(block_length);
                std::iota(output.begin(), output.end(), block_start + 1);
                return output;
            }(),
            by_column,
            ticks,
            map,
            stats
        )
    {}

private:
    DenseCore<solo_, oracle_, Index_, CachedValue_> core;

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        core->fetch_raw(i, buffer);
        return buffer;
    }
};

template<bool oracle_, typename Value_, typename Index_, typename CachedValue_>
struct DenseIndexed : public tatami::DenseExtractor<oracle_, Value_, Index_> {
    DenseIndexed(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& dense_extractor,
        tatami::MaybeOracle<oracle_, Index_> ora,
        tatami::VectorPtr<Index_> indices_ptr,
        bool by_column,
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats& stats) :
        core(
            mat,
            dense_extractor,
            std::move(ora),
            [&]() {
                Rcpp::IntegerVector output(indices_ptr->begin(), indices_ptr->end());
                for (auto& i : output) {
                    ++i;
                }
                return output;
            }(),
            by_column,
            ticks,
            map,
            stats
        )
    {}

private:
    DenseCore<solo_, oracle_, Index_, CachedValue_> core;

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        core->fetch_raw(i, buffer);
        return buffer;
    }
};

}

}

#endif
