#ifndef TATAMI_R_DENSE_EXTRACTOR_HPP
#define TATAMI_R_DENSE_EXTRACTOR_HPP

#include "Rcpp.h"
#include "tatami/tatami.hpp"
#include "tatami_chunked/tatami_chunked.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils.hpp"
#include "parallelize.hpp"
#include "dense_matrix.hpp"

#include <vector>
#include <stdexcept>
#include <type_traits>

namespace tatami_r {

namespace UnknownMatrix_internal {

// GENERAL COMMENTS:
//
// - No need to protect against overflows when incrementing to get to 1-based indexing.
//   This is because the value being incremented is less than the dimension extent, which is known to fit into an Index_.
// - No need to protect against overflows when creating IntegerVectors from dimension extents.
//   We already know that the dimension extent can be safely converted to/from an int, based on checks in the UnknownMatrix constructor.

/********************
 *** Core classes ***
 ********************/

template<bool oracle_, typename Index_> 
class SoloDenseCore {
public:
    SoloDenseCore(
        const Rcpp::RObject& matrix, 
        const Rcpp::Function& dense_extractor,
        bool row,
        tatami::MaybeOracle<oracle_, Index_> oracle,
        Rcpp::IntegerVector non_target_extract, 
        [[maybe_unused]] const std::vector<Index_>& ticks, // provided here for compatibility with the other Dense*Core classes.
        [[maybe_unused]] const std::vector<Index_>& map,
        [[maybe_unused]] const tatami_chunked::SlabCacheStats<Index_>& stats) :
        my_matrix(matrix),
        my_dense_extractor(dense_extractor),
        my_extract_args(2),
        my_row(row),
        my_non_target_length(non_target_extract.size()),
        my_oracle(std::move(oracle))
    {
        my_extract_args[static_cast<int>(row)] = non_target_extract;
    }

private:
    const Rcpp::RObject& my_matrix;
    const Rcpp::Function& my_dense_extractor;
    Rcpp::List my_extract_args;

    bool my_row;
    Index_ my_non_target_length;

    tatami::MaybeOracle<oracle_, Index_> my_oracle;
    typename std::conditional<oracle_, tatami::PredictionIndex, bool>::type my_counter = 0;

public:
    template<typename Value_>
    void fetch_raw(Index_ i, Value_* buffer) {
        if constexpr(oracle_) {
            i = my_oracle->get(my_counter++);
        }

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
        // This involves some Rcpp initializations, so we lock it just in case.
        auto& mexec = executor();
        mexec.run([&]() -> void {
#endif

        my_extract_args[static_cast<int>(!my_row)] = Rcpp::IntegerVector::create(i + 1);
        auto obj = my_dense_extractor(my_matrix, my_extract_args);
        if (my_row) {
            parse_dense_matrix<Index_>(obj, 0, 0, true, buffer, 1, my_non_target_length);
        } else {
            parse_dense_matrix<Index_>(obj, 0, 0, false, buffer, my_non_target_length, 1);
        }

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
        });
#endif
    }
};

template<typename Index_, typename CachedValue_>
class MyopicDenseCore {
public:
    MyopicDenseCore(
        const Rcpp::RObject& matrix, 
        const Rcpp::Function& dense_extractor,
        bool row,
        [[maybe_unused]] tatami::MaybeOracle<false, Index_> oracle, // provided here for compatibility with the other Dense*Core classes.
        Rcpp::IntegerVector non_target_extract, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats<Index_>& stats) :
        my_matrix(matrix),
        my_dense_extractor(dense_extractor),
        my_extract_args(2),
        my_row(row),
        my_non_target_length(non_target_extract.size()),
        my_chunk_ticks(ticks),
        my_chunk_map(map),
        my_factory(stats),
        my_cache(stats.max_slabs_in_cache)
    {
        my_extract_args[static_cast<int>(row)] = non_target_extract;
    }

private:
    const Rcpp::RObject& my_matrix;
    const Rcpp::Function& my_dense_extractor;
    Rcpp::List my_extract_args;

    bool my_row;
    Index_ my_non_target_length;

    const std::vector<Index_>& my_chunk_ticks;
    const std::vector<Index_>& my_chunk_map;

    tatami_chunked::DenseSlabFactory<CachedValue_> my_factory;
    typedef typename decltype(my_factory)::Slab Slab;
    tatami_chunked::LruSlabCache<Index_, Slab> my_cache;

public:
    template<typename Value_>
    void fetch_raw(Index_ i, Value_* buffer) {
        auto chosen = my_chunk_map[i];

        const auto& slab = my_cache.find(
            chosen,
            [&]() -> Slab {
                return my_factory.create();
            },
            [&](Index_ id, Slab& cache) -> void {
#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
                // This involves some Rcpp initializations, so we lock it just in case.
                auto& mexec = executor();
                mexec.run([&]() -> void {
#endif

                auto chunk_start = my_chunk_ticks[id];
                Index_ chunk_len = my_chunk_ticks[id + 1] - chunk_start;
                Rcpp::IntegerVector primary_extract(chunk_len);
                std::iota(primary_extract.begin(), primary_extract.end(), chunk_start + 1);
                my_extract_args[static_cast<int>(!my_row)] = primary_extract;

                auto obj = my_dense_extractor(my_matrix, my_extract_args);
                if (my_row) {
                    parse_dense_matrix<Index_>(obj, 0, 0, true, cache.data, chunk_len, my_non_target_length);
                } else {
                    parse_dense_matrix<Index_>(obj, 0, 0, false, cache.data, my_non_target_length, chunk_len);
                }

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
                });
#endif
            }
        );

        auto shift = sanisizer::product_unsafe<std::size_t>(i - my_chunk_ticks[chosen], my_non_target_length);
        std::copy_n(slab.data + shift, my_non_target_length, buffer);
    }
};

template<typename Index_, typename CachedValue_>
class OracularDenseCore {
public:
    OracularDenseCore(
        const Rcpp::RObject& matrix, 
        const Rcpp::Function& dense_extractor,
        bool row,
        tatami::MaybeOracle<true, Index_> oracle,
        Rcpp::IntegerVector non_target_extract, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats<Index_>& stats) :
        my_matrix(matrix),
        my_dense_extractor(dense_extractor),
        my_extract_args(2),
        my_row(row),
        my_non_target_length(non_target_extract.size()),
        my_chunk_ticks(ticks),
        my_chunk_map(map),
        my_factory(stats),
        my_cache(std::move(oracle), stats.max_slabs_in_cache)
    {
        my_extract_args[static_cast<int>(row)] = non_target_extract;
    }

private:
    const Rcpp::RObject& my_matrix;
    const Rcpp::Function& my_dense_extractor;
    Rcpp::List my_extract_args;

    bool my_row;
    Index_ my_non_target_length;

    const std::vector<Index_>& my_chunk_ticks;
    const std::vector<Index_>& my_chunk_map;

    tatami_chunked::DenseSlabFactory<CachedValue_> my_factory;
    typedef typename decltype(my_factory)::Slab Slab;
    tatami_chunked::OracularSlabCache<Index_, Index_, Slab> my_cache;

public:
    template<typename Value_>
    void fetch_raw(Index_, Value_* buffer) {
        auto res = my_cache.next(
            [&](Index_ i) -> std::pair<Index_, Index_> {
                auto chosen = my_chunk_map[i];
                return std::make_pair(chosen, static_cast<Index_>(i - my_chunk_ticks[chosen]));
            },
            [&]() -> Slab {
                return my_factory.create();
            },
            [&](std::vector<std::pair<Index_, Slab*> >& to_populate) -> void {
                // Sorting them so that the indices are in order.
                auto cmp = [](const std::pair<Index_, Slab*>& left, const std::pair<Index_, Slab*> right) -> bool {
                    return left.first < right.first; 
                };
                if (!std::is_sorted(to_populate.begin(), to_populate.end(), cmp)) {
                    std::sort(to_populate.begin(), to_populate.end(), cmp);
                }

                Index_ total_len = 0;
                for (const auto& p : to_populate) {
                    total_len += my_chunk_ticks[p.first + 1] - my_chunk_ticks[p.first];
                }

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
                // This involves some Rcpp initializations, so we lock it just in case.
                auto& mexec = executor();
                mexec.run([&]() -> void {
#endif

                Rcpp::IntegerVector primary_extract(total_len);
                Index_ current = 0;
                for (const auto& p : to_populate) {
                    Index_ chunk_start = my_chunk_ticks[p.first];
                    Index_ chunk_len = my_chunk_ticks[p.first + 1] - chunk_start;
                    auto start = primary_extract.begin() + current;
                    std::iota(start, start + chunk_len, my_chunk_ticks[p.first] + 1);
                    current += chunk_len;
                }

                my_extract_args[static_cast<int>(!my_row)] = primary_extract;
                auto obj = my_dense_extractor(my_matrix, my_extract_args);

                current = 0;
                for (const auto& p : to_populate) {
                    auto chunk_start = my_chunk_ticks[p.first];
                    Index_ chunk_len = my_chunk_ticks[p.first + 1] - chunk_start;
                    if (my_row) {
                        parse_dense_matrix<Index_>(obj, current, 0, true, p.second->data, chunk_len, my_non_target_length);
                    } else {
                        parse_dense_matrix<Index_>(obj, 0, current, false, p.second->data, my_non_target_length, chunk_len);
                    }
                    current += chunk_len;
                }

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
                });
#endif
            }
        );

        auto shift = sanisizer::product_unsafe<std::size_t>(my_non_target_length, res.second);
        std::copy_n(res.first->data + shift, my_non_target_length, buffer);
    }
};

template<bool solo_, bool oracle_, typename Index_, typename CachedValue_>
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
class DenseFull : public tatami::DenseExtractor<oracle_, Value_, Index_> {
public:
    DenseFull(
        const Rcpp::RObject& matrix, 
        const Rcpp::Function& dense_extractor,
        bool row,
        tatami::MaybeOracle<oracle_, Index_> oracle,
        Index_ non_target_dim,
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats<Index_>& stats) :
        my_core(
            matrix,
            dense_extractor,
            row,
            std::move(oracle),
            [&]{
                Rcpp::IntegerVector output(non_target_dim);
                std::iota(output.begin(), output.end(), 1);
                return output;
            }(),
            ticks,
            map,
            stats
        )
    {}

private:
    DenseCore<solo_, oracle_, Index_, CachedValue_> my_core;

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        my_core.fetch_raw(i, buffer);
        return buffer;
    }
};

template<bool solo_, bool oracle_, typename Value_, typename Index_, typename CachedValue_>
class DenseBlock : public tatami::DenseExtractor<oracle_, Value_, Index_> {
public:
    DenseBlock(
        const Rcpp::RObject& matrix, 
        const Rcpp::Function& dense_extractor,
        bool row,
        tatami::MaybeOracle<oracle_, Index_> oracle,
        Index_ block_start,
        Index_ block_length,
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats<Index_>& stats) :
        my_core(
            matrix,
            dense_extractor,
            row,
            std::move(oracle),
            [&]{
                Rcpp::IntegerVector output(block_length);
                std::iota(output.begin(), output.end(), block_start + 1);
                return output;
            }(),
            ticks,
            map,
            stats
        )
    {}

private:
    DenseCore<solo_, oracle_, Index_, CachedValue_> my_core;

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        my_core.fetch_raw(i, buffer);
        return buffer;
    }
};

template<bool solo_, bool oracle_, typename Value_, typename Index_, typename CachedValue_>
class DenseIndexed : public tatami::DenseExtractor<oracle_, Value_, Index_> {
public:
    DenseIndexed(
        const Rcpp::RObject& matrix, 
        const Rcpp::Function& dense_extractor,
        bool row,
        tatami::MaybeOracle<oracle_, Index_> oracle,
        tatami::VectorPtr<Index_> indices_ptr,
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats<Index_>& stats) :
        my_core(
            matrix,
            dense_extractor,
            row,
            std::move(oracle),
            increment_indices(*indices_ptr),
            ticks,
            map,
            stats
        )
    {}

private:
    DenseCore<solo_, oracle_, Index_, CachedValue_> my_core;

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        my_core.fetch_raw(i, buffer);
        return buffer;
    }
};

}

}

#endif
