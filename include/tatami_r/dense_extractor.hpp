#ifndef TATAMI_R_DENSE_EXTRACTOR_HPP
#define TATAMI_R_DENSE_EXTRACTOR_HPP

#include "Rcpp.h"
#include "tatami/tatami.hpp"
#include "tatami_chunked/tatami_chunked.hpp"
#include "dense_matrix.hpp"

#include <vector>
#include <stdexcept>

namespace tatami_r {

namespace UnknownMatrix_internal {

template<bool oracle_, typename Value_, typename Index_, typename CachedValue_>
struct DenseBase : public tatami::DenseExtractor<oracle_, Value_, Index_> {
    DenseBase(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& dense_extractor,
        tatami::MaybeOracle<oracle_, Index_> ora,
        Rcpp::IntegerVector secondary_extract, 
        bool by_column,
        Index_ max_primary_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        size_t cache_size_in_bytes, 
        bool require_minimum_cache) : 
        mat(mat),
        dense_extractor(dense_extractor),
        extract_args(2),
        by_column(by_column),
        chunk_ticks(ticks),
        chunk_map(map),
        secondary_length(secondary_extract.size()),
        cache(
            max_primary_chunk_length, 
            secondary_length, 
            cache_size_in_bytes / sizeof(CachedValue_), 
            require_minimum_cache, 
            std::move(ora)
        )
    {
        if (cache.num_slabs_in_cache == 0) {
            solo.resize(secondary_length);
        }

        if (by_column) {
            extract_args[0] = secondary_extract;
        } else {
            extract_args[1] = secondary_extract;
        }
    }

    ~DenseBase() = default;

private:
    typedef std::vector<CachedValue_> Slab;

    const Rcpp::RObject& mat;
    const Rcpp::Function& dense_extractor;
    Rcpp::List extract_args;

    bool by_column;
    const std::vector<Index_>& chunk_ticks;
    const std::vector<Index_>& chunk_map;
    size_t secondary_length;

    tatami_chunked::TypicalSlabCacheWorkspace<oracle_, false, Index_, Slab> cache;

    Slab solo;

private:
    std::pair<const Slab*, Index_> fetch_raw(Index_ i) {
        if (cache.num_slabs_in_cache == 0) {
            if constexpr(oracle_) {
                i = cache.cache.next();
            }

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
            // This involves some Rcpp initializations, so we lock it just in case.
            auto& mexec = executor();
            mexec.run([&]() -> void {
#endif

            extract_args[static_cast<int>(by_column)] = Rcpp::IntegerVector::create(i);
            auto obj = dense_extractor(mat, extract_args);
            if (by_column) {
                parse_dense_matrix<false>(obj, solo, 0, 0, secondary_length, 1);
            } else {
                parse_dense_matrix<true>(obj, solo, 0, 0, 1, secondary_length);
            }

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
            });
#endif

            return std::make_pair(&solo, static_cast<Index_>(0));

        } else if constexpr(!oracle_) {
            auto chosen = chunk_map[i];

            const auto& slab = cache.cache.find(
                chosen,
                [&]() -> Slab {
                    return Slab();
                },
                [&](Index_ id, Slab& cache) {
#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
                    // This involves some Rcpp initializations, so we lock it just in case.
                    auto& mexec = executor();
                    mexec.run([&]() -> void {
#endif

                    auto chunk_start = chunk_ticks[id], chunk_end = chunk_ticks[id + 1];
                    size_t chunk_len = chunk_end - chunk_start;
                    Rcpp::IntegerVector primary_extract(chunk_len);
                    std::iota(primary_extract.begin(), primary_extract.end(), chunk_start + 1);
                    extract_args[static_cast<int>(by_column)] = primary_extract;
                    auto obj = dense_extractor(mat, extract_args);
                    if (by_column) {
                        parse_dense_matrix<false>(obj, cache, 0, 0, secondary_length, chunk_len);
                    } else {
                        parse_dense_matrix<true>(obj, cache, 0, 0, chunk_len, secondary_length);
                    }

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
                    });
#endif
                }
            );
            return std::make_pair(&slab, static_cast<Index_>(i - chunk_ticks[chosen]));

        } else {
            return cache.cache.next(
                [&](Index_ i) -> std::pair<Index_, Index_> {
                    auto chosen = chunk_map[i];
                    return std::make_pair(chosen, static_cast<Index_>(i - chunk_ticks[chosen]));
                },
                [&]() -> Slab {
                    return Slab();
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
                            parse_dense_matrix<false>(obj, *p.second, 0, chunk_start, secondary_length, chunk_len);
                        } else {
                            parse_dense_matrix<true>(obj, *p.second, chunk_start, 0, chunk_len, secondary_length);
                        }
                        current += chunk_len;
                    }

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
                    });
#endif
                }
            );
        }
    }

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto res = fetch_raw(i);
        std::copy_n(res.first->data() + this->secondary_length * i, this->secondary_length, buffer);
        return buffer;
    }
};

template<bool oracle_, typename Value_, typename Index_, typename CachedValue_>
struct DenseFull : public DenseBase<oracle_, Value_, Index_, CachedValue_> {
    DenseFull(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& dense_extractor,
        tatami::MaybeOracle<oracle_, Index_> ora,
        Index_ secondary_dim,
        bool by_column,
        Index_ max_primary_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        size_t cache_size_in_bytes, 
        bool require_minimum_cache) :
        DenseBase<oracle_, Value_, Index_, CachedValue_>(
            mat,
            dense_extractor,
            std::move(ora),
            [&]() {
                Rcpp::IntegerVector output(secondary_dim);
                std::iota(output.begin(), output.end(), 1);
                return output;
            }(),
            by_column,
            max_primary_chunk_length,
            ticks,
            map,
            cache_size_in_bytes,
            require_minimum_cache
        )
    {}
};

template<bool oracle_, typename Value_, typename Index_, typename CachedValue_>
struct DenseBlock : public DenseBase<oracle_, Value_, Index_, CachedValue_> {
    DenseBlock(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& dense_extractor,
        tatami::MaybeOracle<oracle_, Index_> ora,
        Index_ block_start,
        Index_ block_length,
        bool by_column,
        Index_ max_primary_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        size_t cache_size_in_bytes, 
        bool require_minimum_cache) :
        DenseBase<oracle_, Value_, Index_, CachedValue_>(
            mat,
            dense_extractor,
            std::move(ora),
            [&]() {
                Rcpp::IntegerVector output(block_length);
                std::iota(output.begin(), output.end(), block_start + 1);
                return output;
            }(),
            by_column,
            max_primary_chunk_length,
            ticks,
            map,
            cache_size_in_bytes,
            require_minimum_cache
        )
    {}
};

template<bool oracle_, typename Value_, typename Index_, typename CachedValue_>
struct DenseIndexed : public DenseBase<oracle_, Value_, Index_, CachedValue_> {
    DenseIndexed(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& dense_extractor,
        tatami::MaybeOracle<oracle_, Index_> ora,
        tatami::VectorPtr<Index_> indices_ptr,
        bool by_column,
        Index_ max_primary_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        size_t cache_size_in_bytes, 
        bool require_minimum_cache) :
        DenseBase<oracle_, Value_, Index_, CachedValue_>(
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
            max_primary_chunk_length,
            ticks,
            map,
            cache_size_in_bytes,
            require_minimum_cache
        )
    {}
};

}

}

#endif
