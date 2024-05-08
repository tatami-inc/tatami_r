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

template<bool oracle_, typename Value_, typename Index_, typename CachedValue_, typename CachedIndex_>
struct SparseBase {
    SparseBase(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& sparse_extractor,
        tatami::MaybeOracle<oracle_, Index_> ora,
        Rcpp::IntegerVector secondary_extract, 
        bool by_column,
        Index_ max_primary_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        size_t cache_size_in_bytes, 
        bool require_minimum_cache,
        bool needs_value,
        bool needs_index) : 
        mat(std::move(mat)),
        sparse_extractor(std::move(sparse_extractor)),
        extract_args(2),
        secondary_indices(std::move(secondary_extract)),
        by_column(by_column),
        chunk_ticks(ticks),
        chunk_map(map),
        max_primary_chunk_length(max_primary_chunk_length),
        secondary_length(secondary_indices.size()),
        needs_value(needs_value),
        needs_index(needs_index),
        cache(
            max_primary_chunk_length, 
            secondary_length, 
            cache_size_in_bytes / std::max(static_cast<size_t>(1), static_cast<size_t>(sizeof(CachedValue_) * needs_value + sizeof(CachedIndex_) * needs_index)),
            require_minimum_cache, 
            std::move(ora)
        )
    {
        if (cache.num_slabs_in_cache == 0) {
            solo = Slab(1, secondary_length, needs_value, needs_index);
        }

        if (by_column) {
            extract_args[0] = secondary_indices;
        } else {
            extract_args[1] = secondary_indices;
        }
    }

    ~SparseBase() = default;

public:
    struct Slab {
        Slab() = default;
        Slab(size_t max_primary_chunk_length, size_t secondary_length, bool needs_value, bool needs_index) {
            if (needs_value) {
                value.reserve(max_primary_chunk_length * secondary_length);
            }
            if (needs_index) {
                index.reserve(max_primary_chunk_length * secondary_length);
            }
            count.reserve(max_primary_chunk_length);
        }

        std::vector<CachedValue_> value;
        std::vector<CachedIndex_> index;
        std::vector<Index_> count;
    };

protected:
    const Rcpp::RObject& mat;
    const Rcpp::Function& sparse_extractor;
    Rcpp::List extract_args;
    Rcpp::IntegerVector secondary_indices;

    bool by_column;
    const std::vector<Index_>& chunk_ticks;
    const std::vector<Index_>& chunk_map;

    size_t max_primary_chunk_length;
    size_t secondary_length;
    bool needs_value;
    bool needs_index;

    std::vector<CachedValue_*> chunk_value_ptrs;
    std::vector<CachedIndex_*> chunk_index_ptrs;
    std::vector<size_t> chunk_counts;

    tatami_chunked::TypicalSlabCacheWorkspace<oracle_, false, Index_, Slab> cache;
    Slab solo;

protected:
    std::pair<const Slab*, Index_> fetch_raw(Index_ i) {
        if (cache.num_slabs_in_cache == 0) {
            if constexpr(oracle_) {
                i = cache.cache.next();
            }

            if (needs_value) {
                chunk_value_ptrs.clear();
                chunk_value_ptrs.push_back(solo.value.data());
            }
            if (needs_index) {
                chunk_index_ptrs.clear();
                chunk_index_ptrs.push_back(solo.index.data());
            }
            solo.count[0] = 0;

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
            // This involves some Rcpp initializations, so we lock it just in case.
            auto& mexec = executor();
            mexec.run([&]() -> void {
#endif

            extract_args[static_cast<int>(by_column)] = Rcpp::IntegerVector::create(i + 1);
            auto obj = sparse_extractor(mat, extract_args);

            if (by_column) {
                parse_sparse_matrix<false>(obj, chunk_value_ptrs, chunk_index_ptrs, solo.count);
            } else {
                parse_sparse_matrix<true>(obj, chunk_value_ptrs, chunk_index_ptrs, solo.count);
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
                    return Slab(max_primary_chunk_length, secondary_length, needs_value, needs_index);
                },
                [&](Index_ id, Slab& cache) {
                    auto chunk_start = chunk_ticks[id], chunk_end = chunk_ticks[id + 1];
                    size_t chunk_len = chunk_end - chunk_start;

                    if (needs_value) {
                        chunk_value_ptrs.clear();
                        auto ptr = cache.value.data();
                        for (size_t i = 0; i < chunk_len; ++i, ptr += secondary_length) {
                            chunk_value_ptrs.push_back(ptr);
                        }
                    }
                    if (needs_index) {
                        chunk_index_ptrs.clear();
                        auto ptr = cache.index.data();
                        for (size_t i = 0; i < chunk_len; ++i, ptr += secondary_length) {
                            chunk_index_ptrs.push_back(ptr);
                        }
                    }
                    cache.count.clear();
                    cache.count.resize(chunk_len);

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
                    // This involves some Rcpp initializations, so we lock it just in case.
                    auto& mexec = executor();
                    mexec.run([&]() -> void {
#endif

                    Rcpp::IntegerVector primary_extract(chunk_len);
                    std::iota(primary_extract.begin(), primary_extract.end(), chunk_start + 1);
                    extract_args[static_cast<int>(by_column)] = primary_extract;
                    auto obj = sparse_extractor(mat, extract_args);

                    if (by_column) {
                        parse_sparse_matrix<false>(obj, chunk_value_ptrs, chunk_index_ptrs, cache.count);
                    } else {
                        parse_sparse_matrix<true>(obj, chunk_value_ptrs, chunk_index_ptrs, cache.count);
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
                    return Slab(max_primary_chunk_length, secondary_length, needs_value, needs_index);
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
                            auto ptr = p.second->value.data();
                            for (Index_ i = 0; i < chunk_len; ++i, ptr += secondary_length) {
                                chunk_value_ptrs.push_back(ptr);
                            }
                        }
                        if (needs_index) {
                            auto ptr = p.second->index.data();
                            for (Index_ i = 0; i < chunk_len; ++i, ptr += secondary_length) {
                                chunk_index_ptrs.push_back(ptr);
                            }
                        }
                    }

                    chunk_counts.clear();
                    chunk_counts.resize(total_len);

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
                    if (by_column) {
                        parse_sparse_matrix<false>(obj, chunk_value_ptrs, chunk_index_ptrs, chunk_counts);
                    } else {
                        parse_sparse_matrix<true>(obj, chunk_value_ptrs, chunk_index_ptrs, chunk_counts);
                    }

                    current = 0;
                    for (const auto& p : to_populate) {
                        Index_ chunk_len = chunk_ticks[p.first + 1] - chunk_ticks[p.first];
                        p.second->count.resize(chunk_len);
                        std::copy_n(chunk_counts.begin() + current, chunk_len, p.second->count.begin());
                        current += chunk_len;
                    }

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
                    });
#endif
                }
            );
        }
    }
};

/******************************
 *** Pure sparse extractors ***
 ******************************/

template<bool oracle_, typename Value_, typename Index_, typename CachedValue_, typename CachedIndex_>
struct SparseFull : public SparseBase<oracle_, Value_, Index_, CachedValue_, CachedIndex_>, public tatami::SparseExtractor<oracle_, Value_, Index_> {
    SparseFull(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& sparse_extractor,
        tatami::MaybeOracle<oracle_, Index_> ora,
        Index_ secondary_dim,
        bool by_column,
        Index_ max_primary_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        size_t cache_size_in_bytes, 
        bool require_minimum_cache,
        bool needs_value,
        bool needs_index) : 
        SparseBase<oracle_, Value_, Index_, CachedValue_, CachedIndex_>(
            std::move(mat),
            std::move(sparse_extractor),
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
            require_minimum_cache,
            needs_value,
            needs_index
        )
    {}

public:
    tatami::SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto res = this->fetch_raw(i);
        const auto& slab = *(res.first);
        Index_ offset = res.second;

        tatami::SparseRange<Value_, Index_> output(slab.count[offset]);
        if (this->needs_value) {
            std::copy_n(slab.value.data() + static_cast<size_t>(offset) * this->secondary_length, output.number, vbuffer); // cast to size_t to avoid overflow.
            output.value = vbuffer;
        }
        if (this->needs_index) {
            std::copy_n(slab.index.data() + static_cast<size_t>(offset) * this->secondary_length, output.number, ibuffer);
            output.index = ibuffer;
        }

        return output;
    }
};

template<bool oracle_, typename Value_, typename Index_, typename CachedValue_, typename CachedIndex_>
struct SparseBlock : public SparseBase<oracle_, Value_, Index_, CachedValue_, CachedIndex_>, public tatami::SparseExtractor<oracle_, Value_, Index_> {
    SparseBlock(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& sparse_extractor,
        tatami::MaybeOracle<oracle_, Index_> ora,
        Index_ block_start,
        Index_ block_length,
        bool by_column,
        Index_ max_primary_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        size_t cache_size_in_bytes, 
        bool require_minimum_cache,
        bool needs_value,
        bool needs_index) : 
        SparseBase<oracle_, Value_, Index_, CachedValue_, CachedIndex_>(
            std::move(mat),
            std::move(sparse_extractor),
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
            require_minimum_cache,
            needs_value,
            needs_index
        ),
        block_start(block_start)
    {}

public:
    tatami::SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto res = this->fetch_raw(i);
        const auto& slab = *(res.first);
        Index_ offset = res.second;

        tatami::SparseRange<Value_, Index_> output(slab.count[offset]);
        if (this->needs_value) {
            std::copy_n(slab.value.data() + static_cast<size_t>(offset) * this->secondary_length, output.number, vbuffer); // cast to size_t to avoid overflow.
            output.value = vbuffer;
        }
        if (this->needs_index) {
            auto iptr = slab.index.data() + static_cast<size_t>(offset) * this->secondary_length;
            for (Index_ i = 0; i < output.number; ++i) {
                ibuffer[i] = static_cast<Index_>(iptr[i]) + block_start;
            }
            output.index = ibuffer;
        }

        return output;
    }

private:
    Index_ block_start;
};

template<bool oracle_, typename Value_, typename Index_, typename CachedValue_, typename CachedIndex_>
struct SparseIndexed : public SparseBase<oracle_, Value_, Index_, CachedValue_, CachedIndex_>, public tatami::SparseExtractor<oracle_, Value_, Index_> {
    SparseIndexed(
        const Rcpp::RObject& mat, 
        const Rcpp::Function& sparse_extractor,
        tatami::MaybeOracle<oracle_, Index_> ora,
        tatami::VectorPtr<Index_> idx_ptr,
        bool by_column,
        Index_ max_primary_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        size_t cache_size_in_bytes, 
        bool require_minimum_cache,
        bool needs_value,
        bool needs_index) : 
        SparseBase<oracle_, Value_, Index_, CachedValue_, CachedIndex_>(
            std::move(mat),
            std::move(sparse_extractor),
            std::move(ora),
            [&]() {
                Rcpp::IntegerVector output(idx_ptr->begin(), idx_ptr->end());
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
            require_minimum_cache,
            needs_value,
            needs_index
        ),
        indices_ptr(std::move(idx_ptr))
    {}

public:
    tatami::SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto res = this->fetch_raw(i);
        const auto& slab = *(res.first);
        Index_ offset = res.second;

        tatami::SparseRange<Value_, Index_> output(slab.count[offset]);
        if (this->needs_value) {
            std::copy_n(slab.value.data() + static_cast<size_t>(offset) * this->secondary_length, output.number, vbuffer); // cast to size_t to avoid overflow.
            output.value = vbuffer;
        }
        if (this->needs_index) {
            auto iptr = slab.index.data() + static_cast<size_t>(offset) * this->secondary_length;
            const auto& indices = *indices_ptr;
            for (Index_ i = 0; i < output.number; ++i) {
                ibuffer[i] = indices[iptr[i]];
            }
            output.index = ibuffer;
        }

        return output;
    }

private:
    tatami::VectorPtr<Index_> indices_ptr;
};

/***********************************
 *** Densified sparse extractors ***
 ***********************************/

template<typename Slab_, typename Value_, typename Index_>
const Value_* densify(const Slab_& slab, Index_ offset, size_t secondary_length, Value_* buffer) {
    size_t shift = static_cast<size_t>(offset) * secondary_length; // cast to size_t to avoid overflow.
    auto vptr = slab.value.data() + shift;
    auto iptr = slab.index.data() + shift;

    std::fill_n(buffer, secondary_length, 0);
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
        Index_ secondary_dim,
        bool by_column,
        Index_ max_primary_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        size_t cache_size_in_bytes, 
        bool require_minimum_cache) :
        SparseBase<oracle_, Value_, Index_, CachedValue_, CachedIndex_>(
            mat,
            sparse_extractor,
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
            require_minimum_cache,
            true,
            true
        )
    {}

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto res = this->fetch_raw(i);
        return densify(*(res.first), res.second, this->secondary_length, buffer);
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
        Index_ max_primary_chunk_length, 
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
            max_primary_chunk_length,
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
        return densify(*(res.first), res.second, this->secondary_length, buffer);
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
        Index_ max_primary_chunk_length, 
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
            max_primary_chunk_length,
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
        return densify(*(res.first), res.second, this->secondary_length, buffer);
    }
};

}

}

#endif
