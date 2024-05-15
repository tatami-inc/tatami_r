#ifndef TATAMI_R_SPARSE_EXTRACTOR_HPP
#define TATAMI_R_SPARSE_EXTRACTOR_HPP

#include "Rcpp.h"
#include "tatami/tatami.hpp"
#include "tatami_chunked/tatami_chunked.hpp"
#include "sparse_matrix.hpp"

#include <vector>
#include <stdexcept>

namespace tatami_r {

namespace UnknownMatrix_internal {

/********************
 *** Core classes ***
 ********************/

template<bool accrow_, bool oracle_, typename Value_, typename Index_>
struct SparseSoloCore {
    SparseSoloCore(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& sparse_extractor,
        tatami::MaybeOracle<oracle_, Index_> ora,
        Rcpp::IntegerVector non_target_extract, 
        [[maybe_unused]] Index_ max_target_chunk_length, // provided here for compatibility with the other Sparse*Core classes.
        [[maybe_unused]] const std::vector<Index_>& ticks,
        [[maybe_unused]] const std::vector<Index_>& map,
        [[maybe_unused]] const tatami_chunked::SlabCacheStats& stats,
        bool needs_value,
        bool needs_index) : 
        mat(mat),
        sparse_extractor(sparse_extractor),
        extract_args(2),
        factory(1, non_target_extract.size(), 1, needs_value, needs_index),
        solo(factory.create())
    {
        extract_args[static_cast<int>(accrow_)] = non_target_extract;
    }

private
    const Rcpp::RObject& mat;
    const Rcpp::Function& sparse_extractor;
    Rcpp::List extract_args;

    tatami_chunked::SparseSlabFactory<CachedValue_, CachedIndex_> factory;
    typedef typename decltype(factory)::Slab Slab;
    Slab solo;

protected:
    std::pair<const Slab*, Index_> fetch_raw(Index_, Value_* vbuffer, Index_* ibuffer) {
        if constexpr(oracle_) {
            i = cache.cache.next();
        }
        solo.number[0] = 0;

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
        // This involves some Rcpp initializations, so we lock it just in case.
        auto& mexec = executor();
        mexec.run([&]() -> void {
#endif

        extract_args[static_cast<int>(!accrow_)] = Rcpp::IntegerVector::create(i + 1);
        auto obj = sparse_extractor(mat, extract_args);
        parse_sparse_matrix<accrow_>(obj, solo.values, solo.indices, solo.number);

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
        });
#endif

        return std::make_pair(&solo, static_cast<Index_>(0));
    }
};

template<bool accrow_, typename Value_, typename Index_, typename CachedValue_, typename CachedIndex_>
struct MyopicSparseCore {
    MyopicSparseCore(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& sparse_extractor,
        [[maybe_unused]] tatami::MaybeOracle<oracle_, Index_> ora, // provided here for compatibility with the other Sparse*Core classes.
        Rcpp::IntegerVector non_target_extract, 
        Index_ max_target_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats& stats,
        bool needs_value,
        bool needs_index) : 
        mat(mat),
        sparse_extractor(sparse_extractor),
        extract_args(2),
        chunk_ticks(ticks),
        chunk_map(map),
        factory(max_target_chunk_length, non_target_extract.size(), stats, needs_value, needs_index),
        cache(stats)
    {
        extract_args[static_cast<int>(accrow_)] = non_target_extract;
    }

private:
    const Rcpp::RObject& mat;
    const Rcpp::Function& sparse_extractor;
    Rcpp::List extract_args;

    const std::vector<Index_>& chunk_ticks;
    const std::vector<Index_>& chunk_map;

    tatami_chunked::SparseSlabFactory<CachedValue_, CachedIndex_> factory;
    typedef typename decltype(factory)::Slab Slab;
    tatami_chunked::LruSlabCache<Index_, Slab> cache;

protected:
    std::pair<const Slab*, Index_> fetch_raw(Index_, Value_* vbuffer, Index_* ibuffer) {
        auto chosen = chunk_map[i];

        const auto& slab = cache.cache.find(
            chosen,
            [&]() -> Slab {
                return factory.create();
            },
            [&](Index_ id, Slab& cache) {
                auto chunk_start = chunk_ticks[id], chunk_end = chunk_ticks[id + 1];
                size_t chunk_len = chunk_end - chunk_start;
                std::fill_n(cache.number, chunk_len, 0);

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
                // This involves some Rcpp initializations, so we lock it just in case.
                auto& mexec = executor();
                mexec.run([&]() -> void {
#endif

                Rcpp::IntegerVector primary_extract(chunk_len);
                std::iota(primary_extract.begin(), primary_extract.end(), chunk_start + 1);
                extract_args[static_cast<int>(!accrow_)] = primary_extract;
                auto obj = sparse_extractor(mat, extract_args);
                parse_sparse_matrix<accrow_>(obj, cache.values, cache.indices, cache.number);

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
                });
#endif
            }
        );

        Index_ offset = i - chunk_ticks[chosen];
        return std::make_pair(&slab, offset);
    }
};

template<bool accrow_, typename Index_, typename CachedValue_, typename CachedIndex_>
struct OracularSparseCore {
    OracularSparseCore(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& sparse_extractor,
        tatami::MaybeOracle<oracle_, Index_> ora,
        Rcpp::IntegerVector non_target_extract, 
        Index_ max_target_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats& stats,
        bool needs_value,
        bool needs_index) : 
        mat(mat),
        sparse_extractor(sparse_extractor),
        extract_args(2),
        chunk_ticks(ticks),
        chunk_map(map),
        factory(max_target_chunk_length, non_target_extract.size(), stats, needs_value, needs_index),
        cache(std::move(ora), stats),
        needs_value(needs_value),
        needs_index(needs_index)
    {
        extract_args[static_cast<int>(accrow_)] = non_target_extract;
    }

private:
    const Rcpp::RObject& mat;
    const Rcpp::Function& sparse_extractor;
    Rcpp::List extract_args;

    const std::vector<Index_>& chunk_ticks;
    const std::vector<Index_>& chunk_map;

    tatami_chunked::SparseSlabFactory<CachedValue_, CachedIndex_> factory;
    typedef typename decltype(factory)::Slab;
    tatami_chunked::OracularSlabCache<Index_, Slab> cache;

    std::vector<CachedValue_*> chunk_value_ptrs;
    std::vector<CachedIndex_*> chunk_index_ptrs;
    std::vector<CachedIndex_> chunk_numbers;

    bool needs_value;
    bool needs_index;

public:
    std::pair<const Slab*, Index_> fetch_raw(Index_) {
        auto res = cache.next(
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

                if (needs_value) {
                    chunk_value_ptrs.clear();
                }
                if (needs_index) {
                    chunk_index_ptrs.clear();
                }

                Index_ total_len = 0;
                for (const auto& p : to_populate) {
                    Index_ chunk_len = chunk_ticks[p.first + 1] - chunk_ticks[p.first];
                    total_len += chunk_len;
                    if (needs_value) {
                        chunk_value_ptrs.insert(chunk_value_ptrs.end(), p.second->values.begin(), p.second->values.end());
                    }
                    if (needs_index) {
                        chunk_index_ptrs.insert(chunk_index_ptrs.end(), p.second->indices.begin(), p.second->indices.end());
                    }
                }

                chunk_numbers.clear();
                chunk_numbers.resize(total_len);

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
                    std::iota(start, start + chunk_len, chunk_start + 1);
                    current += chunk_len;
                }

                extract_args[static_cast<int>(by_column)] = primary_extract;
                auto obj = sparse_extractor(mat, extract_args);
                parse_sparse_matrix<accrow_>(obj, chunk_value_ptrs, chunk_index_ptrs, chunk_numbers.data());

                current = 0;
                for (const auto& p : to_populate) {
                    Index_ chunk_len = chunk_ticks[p.first + 1] - chunk_ticks[p.first];
                    std::copy_n(chunk_numbers.begin() + current, chunk_len, p.second->number);
                    current += chunk_len;
                }

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
                });
#endif
            }
        );
    }
};

template<bool accrow_, bool solo_, bool oracle_, typename Value_, typename Index_, typename CachedValue_, typename CachedIndex_>
using SparseCore = typename std::conditional<solo_,
    SoloSparseCore<accrow_, oracle_, Value_, Index_>,
    typename std::conditional<oracle_,
        OracularSparseCore<accrow_, Index_, CachedValue_, CachedIndex_>,
        MyopicSparseCore<accrow_, Index_, CachedValue_, CachedIndex_>
    >::type
>::type;

/******************************
 *** Pure sparse extractors ***
 ******************************/

template<bool accrow_, bool solo_, bool oracle_, typename Value_, typename Index_, typename CachedValue_, typename CachedIndex_>
struct SparseFull : public tatami::SparseExtractor<oracle_, Value_, Index_> {
    SparseFull(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& sparse_extractor,
        tatami::MaybeOracle<oracle_, Index_> ora,
        Index_ non_target_dim,
        Index_ max_target_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats& stats,
        bool needs_value,
        bool needs_index) : 
        core(
            mat,
            sparse_extractor,
            std::move(ora),
            [&]() {
                Rcpp::IntegerVector output(non_target_dim);
                std::iota(output.begin(), output.end(), 1);
                return output;
            }(),
            max_target_chunk_length,
            ticks,
            map,
            stats,
            needs_value,
            needs_index
        ),
        non_target_dim(non_target_dim),
        needs_value(needs_value),
        needs_index(needs_index)
    {}

private:
    SparseCore<accrow_, solo_, oracle_, Value_, Index_, CachedValue_, CachedIndex_> core;
    Index_ non_target_dim;
    bool needs_value, needs_index;

public:
    tatami::SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto res = core->fetch_raw(i);
        const auto& slab = *(res.first);
        Index_ offset = res.second;

        tatami::SparseRange<Value_, Index_> output(slab.number[offset]);
        if (needs_value) {
            std::copy_n(slab.values[offset], non_target_dim, vbuffer);
            output.value = vbuffer;
        }

        if (needs_index) {
            std::copy_n(slab.indices[offset], non_target_dim, ibuffer);
            output.index = ibuffer;
        }

        return output;
    }
};

template<bool accrow_, bool oracle_, typename Value_, typename Index_, typename CachedValue_, typename CachedIndex_>
struct SparseBlock : public SparseBase<oracle_, Value_, Index_, CachedValue_, CachedIndex_>, public tatami::SparseExtractor<oracle_, Value_, Index_> {
    SparseBlock(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& sparse_extractor,
        tatami::MaybeOracle<oracle_, Index_> ora,
        Index_ block_start,
        Index_ block_length,
        Index_ max_target_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats& stats,
        bool needs_value,
        bool needs_index) : 
        SparseBase<oracle_, Value_, Index_, CachedValue_, CachedIndex_>(
            mat,
            sparse_extractor,
            std::move(ora),
            [&]() {
                Rcpp::IntegerVector output(block_length);
                std::iota(output.begin(), output.end(), block_start + 1);
                return output;
            }(),
            max_target_chunk_length,
            ticks,
            map,
            stats,
            needs_value,
            needs_index
        ),
        block_start(block_start),
        needs_value(needs_value),
        needs_index(needs_index)
    {}

private:
    SparseCore<accrow_, solo_, oracle_, Value_, Index_, CachedValue_, CachedIndex_> core;
    Index_ block_start; 
    bool needs_value, needs_index;

public:
    tatami::SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto res = this->fetch_raw(i);
        const auto& slab = *(res.first);
        Index_ offset = res.second;

        tatami::SparseRange<Value_, Index_> output(slab.count[offset]);
        if (needs_value) {
            std::copy_n(slab.values[offset], output.number, vbuffer); 
            output.value = vbuffer;
        }

        if (needs_index) {
            auto iptr = slab.index[offset];
            for (Index_ i = 0; i < output.number; ++i) {
                ibuffer[i] = static_cast<Index_>(iptr[i]) + block_start;
            }
            output.index = ibuffer;
        }

        return output;
    }
};

template<bool accrow_, bool oracle_, typename Value_, typename Index_, typename CachedValue_, typename CachedIndex_>
struct SparseIndexed : public SparseBase<oracle_, Value_, Index_, CachedValue_, CachedIndex_>, public tatami::SparseExtractor<oracle_, Value_, Index_> {
    SparseIndexed(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& sparse_extractor,
        tatami::MaybeOracle<oracle_, Index_> ora,
        tatami::VectorPtr<Index_> idx_ptr,
        Index_ max_target_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats& stats,
        bool needs_value,
        bool needs_index) : 
        SparseBase<oracle_, Value_, Index_, CachedValue_, CachedIndex_>(
            mat,
            sparse_extractor,
            std::move(ora),
            [&]() {
                Rcpp::IntegerVector output(idx_ptr->begin(), idx_ptr->end());
                for (auto& i : output) {
                    ++i;
                }
                return output;
            }(),
            max_target_chunk_length,
            ticks,
            map,
            stats,
            needs_value,
            needs_index
        ),
        indices_ptr(std::move(idx_ptr))
    {}

private:
    SparseCore<accrow_, solo_, oracle_, Value_, Index_, CachedValue_, CachedIndex_> core;
    tatami::VectorPtr<Index_> indices_ptr;
    bool needs_value, needs_index;

public:
    tatami::SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto res = core->fetch_raw(i);
        const auto& slab = *(res.first);
        Index_ offset = res.second;

        tatami::SparseRange<Value_, Index_> output(slab.count[offset]);
        if (this->needs_value) {
            std::copy_n(slab.values[offset], output.number, vbuffer); 
            output.value = vbuffer;
        }

        if (this->needs_index) {
            auto iptr = slab.indices[offset];
            const auto& indices = *indices_ptr;
            for (CachedIndex_ i = 0; i < output.number; ++i) {
                ibuffer[i] = indices[iptr[i]];
            }
            output.index = ibuffer;
        }

        return output;
    }
};

/***********************************
 *** Densified sparse extractors ***
 ***********************************/

template<typename Slab_, typename Value_, typename Index_>
const Value_* densify(const Slab_& slab, Index_ offset, size_t non_target_length, Value_* buffer) {
    size_t shift = static_cast<size_t>(offset) * non_target_length; // cast to size_t to avoid overflow.
    auto vptr = slab.value.data() + shift;
    auto iptr = slab.index.data() + shift;

    std::fill_n(buffer, non_target_length, 0);
    for (Index_ i = 0, end = slab.count[offset]; i < end; ++i, ++vptr, ++iptr) {
        buffer[*iptr] = *vptr;
    }
    return buffer;
}

template<bool oracle_, typename Value_, typename Index_, typename CachedValue_, typename CachedIndex_>
struct DensifiedSparseFull : public SparseBase<oracle_, Value_, Index_, CachedValue_, CachedIndex_>, public tatami::DenseExtractor<oracle_, Value_, Index_> {
    DensifiedSparseFull(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& sparse_extractor,
        tatami::MaybeOracle<oracle_, Index_> ora,
        Index_ non_target_dim,
        bool by_column,
        Index_ max_target_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        size_t cache_size_in_bytes, 
        bool require_minimum_cache) :
        SparseBase<oracle_, Value_, Index_, CachedValue_, CachedIndex_>(
            mat,
            sparse_extractor,
            std::move(ora),
            [&]() {
                Rcpp::IntegerVector output(non_target_dim);
                std::iota(output.begin(), output.end(), 1);
                return output;
            }(),
            by_column,
            max_target_chunk_length,
            ticks,
            map,
            cache_size_in_bytes,
            require_minimum_cache,
            true,
            true
        )
    {}

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto res = this->fetch_raw(i);
        return densify(*(res.first), res.second, this->non_target_length, buffer);
    }
};

template<bool oracle_, typename Value_, typename Index_, typename CachedValue_, typename CachedIndex_>
struct DensifiedSparseBlock : public SparseBase<oracle_, Value_, Index_, CachedValue_, CachedIndex_>, public tatami::DenseExtractor<oracle_, Value_, Index_> {
    DensifiedSparseBlock(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& sparse_extractor,
        tatami::MaybeOracle<oracle_, Index_> ora,
        Index_ block_start,
        Index_ block_length,
        bool by_column,
        Index_ max_target_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        size_t cache_size_in_bytes, 
        bool require_minimum_cache) :
        SparseBase<oracle_, Value_, Index_, CachedValue_, CachedIndex_>(
            mat,
            sparse_extractor,
            std::move(ora),
            [&]() {
                Rcpp::IntegerVector output(block_length);
                std::iota(output.begin(), output.end(), block_start + 1);
                return output;
            }(),
            by_column,
            max_target_chunk_length,
            ticks,
            map,
            cache_size_in_bytes,
            require_minimum_cache,
            true,
            true
        )
    {}

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto res = this->fetch_raw(i);
        return densify(*(res.first), res.second, this->non_target_length, buffer);
    }
};

template<bool oracle_, typename Value_, typename Index_, typename CachedValue_, typename CachedIndex_>
struct DensifiedSparseIndexed : public SparseBase<oracle_, Value_, Index_, CachedValue_, CachedIndex_>, public tatami::DenseExtractor<oracle_, Value_, Index_> {
    DensifiedSparseIndexed(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& sparse_extractor,
        tatami::MaybeOracle<oracle_, Index_> ora,
        tatami::VectorPtr<Index_> idx_ptr,
        bool by_column,
        Index_ max_target_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        size_t cache_size_in_bytes, 
        bool require_minimum_cache) :
        SparseBase<oracle_, Value_, Index_, CachedValue_, CachedIndex_>(
            mat,
            sparse_extractor,
            std::move(ora),
            [&]() {
                Rcpp::IntegerVector output(idx_ptr->begin(), idx_ptr->end());
                for (auto& i : output) {
                    ++i;
                }
                return output;
            }(),
            by_column,
            max_target_chunk_length,
            ticks,
            map,
            cache_size_in_bytes,
            require_minimum_cache,
            true,
            true
        )
    {}

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto res = this->fetch_raw(i);
        return densify(*(res.first), res.second, this->non_target_length, buffer);
    }
};

}

}

#endif
