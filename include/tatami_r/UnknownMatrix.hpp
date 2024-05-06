#ifndef TATAMI_R_UNKNOWNMATRIX_HPP
#define TATAMI_R_UNKNOWNMATRIX_HPP

#include "Rcpp.h"
#include "tatami/tatami.hpp"
#include "tatami_chunked/tatami_chunked.hpp"
#include "SimpleMatrix.hpp"
#include "SVT_SparseMatrix.hpp"
#include "COO_SparseMatrix.hpp"

#include <vector>
#include <memory>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include "parallelize.hpp"

namespace tatami_r {

/**
 * @cond
 */
namespace UnknownMatrix_internal {

/*****************
 ***** Dense *****
 *****************/

template<bool oracle_, typename Index_, typename CachedValue_>
struct DenseBase {
    DenseBase(
        Rcpp::RObject mat, 
        Rcpp::Function dense_extractor,
        Index_ primary_chunk_length, 
        Rcpp::IntegerVector secondary_extract, 
        bool by_column,
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        size_t cache_size_in_bytes, 
        bool require_minimum_cache,
        tatami::MaybeOracle<oracle_, Index_> ora) : 
        mat(std::move(mat)),
        dense_extractor(std::move(dense_extractor)),
        extract_args(2),
        by_column(by_column),
        chunk_ticks(ticks),
        chunk_map(map),
        cache(primary_chunk_length, secondary_extract.size(), cache_size_in_elements / sizeof(CachedValue_), require_minimum_cache, std::move(ora))
    {
        if (cache.num_slabs_in_cache == 0) {
            solo.resize(secondary_extract.size());
        }

        if (by_column_) {
            extract_args[0] = secondary_extract;
        } else {
            extract_args[1] = secondary_extract;
        }
    }

    std::pair<const Slab*, Index_> fetch(Index_ i) {
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
            parse_simple_matrix(obj, solo, !by_column);

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
                    Rcpp::IntegerVector primary_extract(chunk_end - chunk_start);
                    std::iota(primary_extract.begin(), primary_extract.end(), chunk_start + 1);
                    extract_args[static_cast<int>(by_column)] = primary_extract;
                    auto obj = dense_extractor(mat, extract_args);
                    parse_simple_matrix(obj, cache, !by_column);

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

                    Rcpp::IntegerVector primary_extract(len);
                    Index_ current = 0;
                    for (const auto& p : to_populate) {
                        Index_ to_fill = chunk_ticks[p.first + 1] - chunk_ticks[p.first];
                        auto start = primary_extract.begin() + current;
                        std::iota(start, start + to_fill, chunk_ticks[p.first] + 1);
                        current += to_fill;
                    }

                    extract_args[static_cast<int>(by_column)] = primary_extract;
                    auto obj = dense_extractor(mat, extract_args);

                    current = 0;
                    for (const auto& p : to_populate) {
                        Index_ to_fill = chunk_ticks[p.first + 1] - chunk_ticks[p.first];
                        parse_simple_matrix(obj, cache, transpose, current, to_fill);
                        current += to_fill;
                    }

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
                    });
#endif
                }
            );
        }
    }

public:
    typedef std::vector<CachedValue_> Slab;

protected:
    Rcpp::RObject mat;
    Rcpp::Function dense_extractor;
    Rcpp::List extract_args;

    bool by_column;
    const std::vector<Index_>& chunk_ticks;
    const std::vector<Index_>& chunk_map;

    tatami_chunked::TypicalCacheWorkspace<oracle_, false, Index_, Slab> cache;

    Slab solo;
};

template<bool oracle_, typename Value_, typename Index_, typename CachedValue_>
struct DenseFull : public DenseBase<oracle_, Index_, CachedValue_> {
    DenseFull(
        Rcpp::RObject mat, 
        Rcpp::Function dense_extractor,
        Index_ primary_chunk_length, 
        Index_ secondary_dim,
        bool by_column,
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        size_t cache_size_in_bytes, 
        bool require_minimum_cache,
        tatami::MaybeOracle<oracle_, Index_> ora) : 
        DenseBase<oracle_, Value_, Index_, CachedValue_>(
            std::move(mat),
            std::move(dense_extractor),
            primary_chunk_length,
            [&]() {
                Rcpp::IntegerVector output(secondary_dim);
                std::iota(output.begin(), output.end(), 1);
                return output;
            }(),
            by_column,
            ticks,
            map,
            cache_size_in_bytes,
            require_minimum_cache
        ),
        full_extent(secondary_dim);
    {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto res = fetch_raw(i);
        std::copy(res.first.data() + full_extent * i, full_extent, buffer);
        return buffer;
    }

private:
    size_t full_extent;
};

template<bool oracle_, typename Value_, typename Index_, typename CachedValue_>
struct DenseBlock : public DenseBase<oracle_, Index_, CachedValue_> {
    DenseBlock(
        Rcpp::RObject mat, 
        Rcpp::Function dense_extractor,
        Index_ primary_chunk_length, 
        Index_ block_start,
        Index_ block_length,
        bool by_column,
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        size_t cache_size_in_bytes, 
        bool require_minimum_cache,
        tatami::MaybeOracle<oracle_, Index_> ora) : 
        DenseBase<oracle_, Value_, Index_, CachedValue_>(
            std::move(mat),
            std::move(dense_extractor),
            primary_chunk_length,
            [&]() {
                Rcpp::IntegerVector output(block_length);
                std::iota(output.begin(), output.end(), block_start + 1);
                return output;
            }(),
            by_column,
            ticks,
            map,
            cache_size_in_bytes,
            require_minimum_cache
        ),
        block_length(secondary_dim);
    {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto res = fetch_raw(i);
        std::copy(res.first.data() + block_length * i, block_length, buffer);
        return buffer;
    }

private:
    size_t block_length;
};

template<bool oracle_, typename Value_, typename Index_, typename CachedValue_>
struct DenseIndexed : public DenseBase<oracle_, Index_, CachedValue_> {
    DenseIndexed(
        Rcpp::RObject mat, 
        Rcpp::Function dense_extractor,
        Index_ primary_chunk_length, 
        const std::vector<Index_>& indices,
        bool by_column,
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        size_t cache_size_in_bytes, 
        bool require_minimum_cache,
        tatami::MaybeOracle<oracle_, Index_> ora) : 
        DenseBase<oracle_, Value_, Index_, CachedValue_>(
            std::move(mat),
            std::move(dense_extractor),
            primary_chunk_length,
            [&]() {
                Rcpp::IntegerVector output(indices.begin(), indices.end());
                for (auto& i : output) {
                    ++i;
                }
                return output;
            }(),
            by_column,
            ticks,
            map,
            cache_size_in_bytes,
            require_minimum_cache
        ),
        num_indices(indices.size())
    {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto res = fetch_raw(i);
        std::copy(res.first.data() + num_indices * i, num_indices, buffer);
        return buffer;
    }

private:
    size_t num_indices;
};

/******************
 ***** Sparse *****
 ******************/

}
/**
 * @endcond
 */


/**
 * @brief Unknown matrix-like object in R.
 *
 * @tparam Value_ Numeric type of data value for the interface.
 * @tparam Index_ Integer type for the row/column indices, for the interface.
 *
 * Pull data out of an unknown matrix-like object by calling methods from the [**DelayedArray**](https://bioconductor.org/packages/DelayedArray) package via **Rcpp**.
 * This effectively extends **tatami** to work with any abstract numeric matrix that might be consumed by an R function.
 * Note that this class should only be constructed in a serial context, and some additional effort is required to deal with parallelization of its methods; see `executor()` for more details.
 */
template<typename Value_, typename Index_>
class UnknownMatrix : public tatami::Matrix<Value_, Index_> {
public:
    /**
     * @param seed A matrix-like R object.
     * @param cache Size of the cache, in bytes.
     * If -1, this is determined from `DelayedArray::getAutoBlockSize()`.
     *
     * This constructor should only be called in a serial context, as the (default) construction of **Rcpp** objects may call the R API.
     */
    UnknownMatrix(Rcpp::RObject seed, size_t cache = -1) : 
        original_seed(seed), 
        delayed_env(Rcpp::Environment::namespace_env("DelayedArray")),
        sparse_env(Rcpp::Environment::namespace_env("SparseArray")),
        dense_extractor(delayed_env["extract_array"]),
        sparse_extractor(sparse_env["extract_sparse_array"])
    {
        // We assume the constructor only occurs on the main thread, so we
        // won't bother locking things up. I'm also not sure that the
        // operations in the initialization list are thread-safe.

        {
            auto base = Rcpp::Environment::base_env();
            Rcpp::Function fun = base["dim"];
            Rcpp::RObject output = fun(seed);
            if (output.sexp_type() != INTSXP) {
                auto ctype = get_class_name(original_seed);
                throw std::runtime_error("'dim(<" + ctype + ">)' should return an integer vector");
            }
            Rcpp::IntegerVector dims(output);
            if (dims.size() != 2 || dims[0] < 0 || dims[1] < 0) {
                auto ctype = get_class_name(original_seed);
                throw std::runtime_error("'dim(<" + ctype + ">)' should contain two non-negative integers");
            }

            internal_nrow = dims[0];
            internal_ncol = dims[1];
        }

        {
            Rcpp::Function fun = delayed_env["is_sparse"];
            Rcpp::LogicalVector is_sparse = fun(seed);
            if (is_sparse.size() != 1) {
                auto ctype = get_class_name(original_seed);
                throw std::runtime_error("'is_sparse(<" + ctype + ">)' should return a logical vector of length 1");
            }
            internal_sparse = (is_sparse[0] != 0);
        }

        {
            row_chunk_map.resize(internal_nrow);
            col_chunk_map.resize(internal_ncol);

            Rcpp::Function fun = delayed_env["chunkGrid"];
            Rcpp::RObject output = fun(seed);

            if (output == R_NilValue) {
                row_max_chunk_size = 1;
                col_max_chunk_size = 1;
                std::iota(row_chunk_map.begin(), row_chunk_map.end(), static_cast<Index_>(0));
                std::iota(col_chunk_map.begin(), col_chunk_map.end(), static_cast<Index_>(0));
                row_chunk_ticks.resize(internal_nrow + 1);
                std::iota(row_chunk_ticks.begin(), row_chunk_ticks.end(), static_cast<Index_>(0));
                col_chunk_ticks.resize(internal_col + 1);
                std::iota(col_chunk_ticks.begin(), col_chunk_ticks.end(), static_cast<Index_>(0));

            } else {
                auto grid_cls = get_class_name(output);

                if (grid_cls == "RegularArrayGrid") {
                    Rcpp::IntegerVector spacings(Rcpp::RObject(object.slot("spacings")));
                    if (spacings.size() != 2) {
                        throw std::runtime_error("'chunkGrid(<" + ctype + ">)@spacings' should be an integer vector of length 2 with non-negative values");
                    }

                    auto populate = [](Index_ extent, Index_ spacing, std::vector<Index_>& map, std::vector<Index_>& chunks) {
                        ticks.reserve((extent / spacing) + (extent % spacing > 0) + 1);
                        Index_ start = 0;
                        ticks.push_back(start);
                        while (start != extent) {
                            auto to_fill = std::min(spacing, extent - start);
                            std::fill_n(map.begin() + start, to_fill, chunks.size());
                            start += to_fill;
                            ticks.push_back(start);
                        }
                    };

                    row_max_chunk_size = spacings[0];
                    if (row_max_chunk_size > 0) {
                        populate(internal_nrow, row_max_chunk_size, row_chunk_map, row_chunk_ticks);
                    }
                    col_max_chunk_size = chunk_sizes[1];
                    if (col_max_chunk_size > 0) {
                        populate(internal_ncol, col_max_chunk_size, col_chunk_map, col_chunk_ticks);
                    }

                } else if (grid_cls == "ArbitraryArrayGrid") {
                    Rcpp::List ticks(Rcpp::RObject(object.slot("tickmarks")));
                    if (ticks.size() != 2) {
                        throw std::runtime_error("'chunkGrid(<" + ctype + ">)@tickmarks' should return a list of length 2");
                    }

                    auto populate = [](Index_ extent, const Rcpp::IntegerVector& ticks, std::vector<Index_>& map, std::vector<Index_>& new_ticks, Index_& chunk_size) {
                        if (ticks.size() == 0 || ticks.back() != static_cast<int>(extent)) {
                            throw std::runtime_error("invalid ticks in 'chunkGrid(<" + ctype + ">)@tickmarks");
                        }
                        std::copy(ticks.begin(), ticks.end(), new_ticks.begin());

                        max_chunk_size = 0;
                        chunks.reserve(ticks.size());
                        int start = 0;

                        for (auto t : ticks) {
                            if (t < start) {
                                throw std::runtime_error("invalid ticks in 'chunkGrid(<" + ctype + ">)@tickmarks'");
                            }
                            auto to_fill = t - start;
                            if (to_fill > chunk_size) {
                                max_chunk_size = to_fill;
                            }
                            std::fill_n(row_chunk_map.begin() + start, to_fill, row_chunks.size());
                            start = t;
                        }
                    };

                    Rcpp::IntegerVector first(ticks[0]);
                    populate(internal_nrow, first, row_chunk_map, row_chunks, row_max_chunk_size);
                    Rcpp::IntegerVector second(ticks[1]);
                    populate(internal_ncol, second, col_chunk_map, col_chunks, col_max_chunk_size);

                } else {
                    throw std::runtime_error("instance of unknown class '" + grid_cls + "' returned by 'chunkGrid(<" + ctype + ">)");
                }
            }
        }

        cache_size_in_bytes = cache;
        if (cache_size_in_bytes == static_cast<size_t>(-1)) {
            Rcpp::Function fun = delayed_env["getAutoBlockSize"];
            Rcpp::NumericVector output = fun();
            if (output.size() != 1 || output[0] < 0) {
                throw std::runtime_error("'getAutoBlockSize()' should return a non-negative number of bytes");
            }
            cache_size_in_bytes = output[0];
        }

        auto chunks_per_row = col_chunks.size();
        auto chunks_per_col = row_chunks.size();
        internal_prefer_rows = chunks_per_row <= chunks_per_col;
    }

private:
    Index_ internal_nrow, internal_ncol;
    bool internal_sparse, internal_prefer_rows;

    std::vector<Index_> row_chunk_map, col_chunk_map;
    std::vector<Index_> row_chunk_ticks, col_chunk_ticks;

    // To decide how many chunks to store in the cache, we pretend the largest
    // chunk is a good representative. This is a bit suboptimal for irregular
    // chunks but the LruSlabCache class doesn't have a good way of dealing
    // with this right now. The fundamental problem is that variable slabs will
    // either (i) all reach the maximum allocation eventually, if slabs are
    // reused, or (ii) require lots of allocations, if slabs are not reused, or
    // (iii) require manual defragmentation, if slabs are reused in a manner
    // that avoids inflation to the maximum allocation.
    Index_ row_max_chunk_size, col_max_chunk_size;

    size_t cache_size_in_bytes;

    Rcpp::RObject original_seed;
    Rcpp::Environment delayed_env, sparse_env;
    Rcpp::Function dense_extractor, sparse_extractor;

public:
    Index_ nrow() const {
        return internal_nrow;
    }

    Index_ ncol() const {
        return internal_ncol;
    }

    bool sparse() const {
        return internal_sparse;
    }

    double sparse_proportion() const {
        return static_cast<double>(internal_sparse);
    }

    bool prefer_rows() const {
        return internal_prefer_rows;
    }

    double prefer_rows_proportion() const {
        return static_cast<double>(internal_prefer_rows);
    }

    bool uses_oracle(bool) const {
        return true;
    }

private:
    template<bool sparse_>
    struct Workspace {
        Workspace(Index_ full) {
            secondary_indices = R_NilValue;
            secondary_len = full;
        }

        Workspace(Index_ start, Index_ length) {
            secondary_indices = create_consecutive_indices(start, length);
            secondary_len = length;
        }

        Workspace(const std::vector<Index_>& indices) {
            Rcpp::IntegerVector temp(indices.begin(), indices.end());
            for (auto& x : temp) { ++x; } // 1-based
            secondary_indices = temp;
            secondary_len = indices.size();
        }

    public:
        Rcpp::RObject secondary_indices;
        Index_ secondary_len;

        std::shared_ptr<tatami::Matrix<Value_, Index_> > buffer;
        std::shared_ptr<tatami::Extractor<tatami::DimensionSelectionType::FULL, sparse_, Value_, Index_> > bufextractor;
        Rcpp::RObject contents;

        // We define the number of chunks in the extracted block,
        // as well as the size of the extracted block itself.
        size_t max_block_size_in_chunks;
        size_t max_block_size_in_elements;

        // No oracle, we just go for blocks.
        Index_ primary_block_start, primary_block_len;

        // With an oracle, we try to make some better decisions.
        tatami::OracleStream<Index_> prediction_stream;
        std::unordered_map<Index_, Index_> mapping;
        std::unordered_set<Index_> used_chunks;
        std::vector<Index_> unique_predictions;
        size_t predictions_used = 0, predictions_made = 0;
    };

private:
    static std::pair<Index_, Index_> round_indices(Index_ i, Index_ interval, Index_ max) {
        Index_ new_first = (i / interval) * interval;
        Index_ new_last = std::min(max, new_first + interval);
        return std::make_pair(new_first, new_last - new_first);
    }

    static Rcpp::IntegerVector create_consecutive_indices(Index_ start, Index_ length) {
        Rcpp::IntegerVector out(length);
        std::iota(out.begin(), out.end(), start + 1); // 1-based
        return out;
    }

    template<bool byrow_, bool sparse_>
    Rcpp::List create_rounded_indices(Index_ i, Workspace<sparse_>* work) const {
        Rcpp::List indices(2);
        if constexpr(byrow_) {
            auto row_rounded = round_indices(i, work->max_block_size_in_elements, internal_nrow);
            work->primary_block_start = row_rounded.first;
            work->primary_block_len = row_rounded.second;
            indices[0] = create_consecutive_indices(row_rounded.first, row_rounded.second);
            indices[1] = work->secondary_indices;
        } else {
            auto col_rounded = round_indices(i, work->max_block_size_in_elements, internal_ncol);
            work->primary_block_start = col_rounded.first;
            work->primary_block_len = col_rounded.second;
            indices[0] = work->secondary_indices;
            indices[1] = create_consecutive_indices(col_rounded.first, col_rounded.second);
        }
        return indices;
    }

    template<bool byrow_>
    Index_ chunk_dim() const {
        if constexpr(byrow_) {
            return chunk_nrow;
        } else {
            return chunk_ncol;
        }
    }

    template<bool byrow_, bool sparse_>
    Rcpp::List create_next_indices(Workspace<sparse_>* work) const {
        auto& up = work->unique_predictions;
        auto& map = work->mapping;
        auto& used = work->used_chunks;
        up.clear();
        map.clear();
        used.clear();

        auto& stream = work->prediction_stream;
        work->predictions_used = 0;
        auto& made = work->predictions_made;
        made = 0;
        size_t max_predictions = work->max_block_size_in_elements * 10;
        auto cdim = chunk_dim<byrow_>();

        for (size_t i = 0; i < max_predictions; ++i) {
            Index_ current;
            if (!stream.next(current)) {
                break;
            }

            auto it = map.find(current);
            if (it == map.end()) {
                // If the chunking is non-trivial, we quit when the cached
                // number of chunks is reached. This is equivalent to only
                // loading an element from a chunk if we're sure that we can
                // cache the entire chunk (even if we don't actually load the
                // entire chunk), so as to avoid running out of space halfway
                // through the chunk and requiring another read of the same
                // chunk in the next set of predictions.
                if (cdim > 1) {
                    Index_ chunk = current / cdim;
                    auto uIt = used.find(chunk);
                    if (uIt == used.end()) {
                        if (used.size() == work->max_block_size_in_chunks) {
                            stream.back();
                            break;
                        }
                        used.insert(chunk);
                    }

                } else if (map.size() == work->max_block_size_in_elements) {
                    stream.back();
                    break;                    
                }

                map[current] = up.size();
                up.push_back(current);
            }

            ++made;
        }

        // Creating sorted, 1-based indices to use in the DelayedArray extractors.
        if (!std::is_sorted(up.begin(), up.end())) {
            std::sort(up.begin(), up.end());
            Index_ counter = 0;
            for (auto x : up) {
                map[x] = counter;
                ++counter;
            }
        }

        Rcpp::IntegerVector primary_indices(up.begin(), up.end());
        for (auto& x : primary_indices) { 
            ++x;
        }
        work->primary_block_len = primary_indices.size();

        Rcpp::List indices(2);
        if constexpr(byrow_) {
            indices[0] = std::move(primary_indices);
            indices[1] = work->secondary_indices;
        } else {
            indices[0] = work->secondary_indices;
            indices[1] = std::move(primary_indices);
        }
        return indices;
    }

    template<bool byrow_, bool sparse_, bool sparse_err = sparse_>
    void check_buffered_dims(const tatami::Matrix<Value_, Index_>* parsed, const Workspace<sparse_>* work) const {
        Index_ parsed_primary = (byrow_ ? parsed->nrow() : parsed->ncol());
        Index_ parsed_secondary = (byrow_ ? parsed->ncol() : parsed->nrow());

        if (parsed_primary != work->primary_block_len || parsed_secondary != work->secondary_len) {
            auto ctype = get_class_name(original_seed);
            throw std::runtime_error("'" + 
                (sparse_err ? std::string("extract_sparse_array") : std::string("extract_array")) + 
                "(<" + ctype + ">)' returns incorrect dimensions");
        }
    }

    template<bool sparse_>
    static bool needs_reset(Index_ i, const Workspace<sparse_>* work) {
        if (work->buffer != nullptr) {
            if (i >= work->primary_block_start && i < work->primary_block_start + work->primary_block_len) {
                return false;
            }
        }
        return true;
    }

private:
    template<bool byrow_, class Function_>
    void run_dense_extractor(Index_ i, const tatami::Options& options, Workspace<false>* work, Function_ index_creator) const {
#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
        // This involves some Rcpp initializations, so we lock it just in case.
        auto& mexec = executor();
        mexec.run([&]() -> void {
#endif

        auto indices = index_creator(i, work);
        Rcpp::RObject val0 = dense_extractor(original_seed, indices);
        auto parsed = parse_simple_matrix<Value_, Index_>(val0);
        work->contents = std::move(parsed.contents);
        work->buffer = std::move(parsed.matrix);

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
        });
#endif

        check_buffered_dims<byrow_, false>(work->buffer.get(), work);
        work->bufextractor = tatami::new_extractor<byrow_, false>(work->buffer.get(), options);
    }

    template<bool byrow_>
    const Value_* run_dense_extractor(Index_ i, Value_* buffer, const tatami::Options& options, Workspace<false>* work) const {
        if (work->prediction_stream.active()) {
            if (work->predictions_used == work->predictions_made) {
                run_dense_extractor<byrow_>(i, options, work, [&](Index_, Workspace<false>* work2) -> Rcpp::List { 
                    return this->create_next_indices<byrow_>(work2);
                });
            }
            i = work->mapping.find(i)->second;
            ++work->predictions_used;

        } else {
            if (needs_reset(i, work)) {
                run_dense_extractor<byrow_>(i, options, work, [&](Index_ i2, Workspace<false>* work2) -> Rcpp::List { 
                    return this->create_rounded_indices<byrow_>(i2, work2);
                });
            }
            i -= work->primary_block_start;
        }

        // Forcing a copy to avoid a possible reference to a transient buffer inside 'work->buffer'.
        return work->bufextractor->fetch_copy(i, buffer);
    }

    template<bool byrow_, class Function_>
    void run_sparse_extractor(Index_ i, const tatami::Options& options, Workspace<true>* work, Function_ index_creator) const {
#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
        // This involves some Rcpp initializations, so we lock it just in case.
        auto& mexec = executor();
        mexec.run([&]() -> void {
#endif

        auto indices = index_creator(i, work);

        if (internal_sparse) {
            auto val0 = sparse_extractor(original_seed, indices);

            auto ctype = get_class_name(val0);
            Parsed<Value_, Index_> parsed;
            if (ctype == "SVT_SparseMatrix") {
                parsed = parse_SVT_SparseMatrix<Value_, Index_>(val0);
            } else if (ctype == "COO_SparseMatrix") {
                parsed = parse_COO_SparseMatrix<Value_, Index_>(val0, byrow_);
            } else {
                throw std::runtime_error("unknown class '" + ctype + "' returned from 'extract_sparse_array()'");
            }
            check_buffered_dims<byrow_, true, true>(parsed.matrix.get(), work);

            work->buffer = std::move(parsed.matrix);
            work->contents = std::move(parsed.contents);
        } else {
            auto val0 = dense_extractor(original_seed, indices);
            auto parsed = parse_simple_matrix<Value_, Index_>(val0);
            check_buffered_dims<byrow_, true, false>(parsed.matrix.get(), work);

            work->buffer = std::move(parsed.matrix);
            work->contents = std::move(parsed.contents);
        }

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
        });
#endif

        work->bufextractor = tatami::new_extractor<byrow_, true>(work->buffer.get(), options);
    }

    template<bool byrow_>
    tatami::SparseRange<Value_, Index_> run_sparse_extractor(Index_ i, Value_* vbuffer, Index_* ibuffer, const tatami::Options& options, Workspace<true>* work) const {
        if (work->prediction_stream.active()) {
            if (work->predictions_used == work->predictions_made) {
                run_sparse_extractor<byrow_>(i, options, work, [&](Index_, Workspace<true>* work2) -> Rcpp::List { 
                    return this->create_next_indices<byrow_>(work2);
                });
            }
            i = work->mapping.find(i)->second;
            ++work->predictions_used;

        } else {
            if (needs_reset(i, work)) {
                run_sparse_extractor<byrow_>(i, options, work, [&](Index_ i2, Workspace<true>* work2) -> Rcpp::List { 
                    return this->create_rounded_indices<byrow_>(i2, work2);
                });
            }
            i -= work->primary_block_start;
        }

        // Forcing a copy to avoid a possible reference to a transient buffer inside 'work->buffer'.
        return work->bufextractor->fetch_copy(i, vbuffer, ibuffer);
    }

private:
    template<bool byrow_, tatami::DimensionSelectionType selection_, bool sparse_>
    struct UnknownExtractor : public tatami::Extractor<selection_, sparse_, Value_, Index_> {
        template<typename ... Args_>
        static auto setup_workspace(Args_&&... args) {
            Workspace<sparse_>* tmp;

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
            // This involves some Rcpp initializations, so we lock it just in case.
            auto& mexec = executor();
            mexec.run([&]() -> void {
#endif
            tmp = new Workspace<sparse_>(std::forward<Args_>(args)...);
#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
            });
#endif

            return tmp;
        }

        static void define_block_size(const UnknownMatrix<Value_, Index_>* parent, Index_ len, Workspace<sparse_>& work) {
            double cache_elements = static_cast<double>(parent->cache_size) / (static_cast<double>(len) * static_cast<double>(sizeof(Value_)));
            auto chunk_dim = parent->template chunk_dim<byrow_>();
            work.max_block_size_in_chunks = std::max(1.0, std::floor(cache_elements / chunk_dim));
            work.max_block_size_in_elements = work.max_block_size_in_chunks * chunk_dim;
        }

        UnknownExtractor(const UnknownMatrix<Value_, Index_>* p) : parent(p) { 
            if constexpr(selection_ == tatami::DimensionSelectionType::FULL) {
                this->full_length = byrow_ ? p->internal_ncol : p->internal_nrow;
                work.reset(setup_workspace(this->full_length));
                define_block_size(parent, this->full_length, *work);
            }
        }

        UnknownExtractor(const UnknownMatrix<Value_, Index_>* p, Index_ start, Index_ length) : parent(p) {
            if constexpr(selection_ == tatami::DimensionSelectionType::BLOCK) {
                this->block_start = start;
                this->block_length = length;
                work.reset(setup_workspace(start, length));
                define_block_size(parent, this->block_length, *work);
            }
        }

        UnknownExtractor(const UnknownMatrix<Value_, Index_>* p, std::vector<Index_> idx) : parent(p) {
            if constexpr(selection_ == tatami::DimensionSelectionType::INDEX) {
                indices = std::move(idx);
                this->index_length = indices.size();
                work.reset(setup_workspace(indices));
                define_block_size(parent, this->index_length, *work);
            }
        }

        const Index_* index_start() const {
            if constexpr(selection_ == tatami::DimensionSelectionType::INDEX) {
                return indices.data();
            } else {
                return NULL;
            }
        }

    protected:
        const UnknownMatrix<Value_, Index_>* parent;
        std::unique_ptr<Workspace<sparse_> > work;
        typename std::conditional<selection_ == tatami::DimensionSelectionType::INDEX, std::vector<Index_>, bool>::type indices;
    };

private:
    template<bool byrow_, tatami::DimensionSelectionType selection_>
    struct DenseUnknownExtractor : public UnknownExtractor<byrow_, selection_, false> {
        template<typename ... Args_>
        DenseUnknownExtractor(const UnknownMatrix<Value_, Index_>* p, tatami::Options opt, Args_&&... args) : 
            UnknownExtractor<byrow_, selection_, false>(p, std::forward<Args_>(args)...), options(std::move(opt)) {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            return this->parent->template run_dense_extractor<byrow_>(i, buffer, options, this->work.get());
        }

    private:
        tatami::Options options;
    };

private:
    template<bool byrow_, tatami::DimensionSelectionType selection_>
    struct SparseUnknownExtractor : public UnknownExtractor<byrow_, selection_, true> {
        template<typename ... Args_>
        SparseUnknownExtractor(const UnknownMatrix<Value_, Index_>* p, tatami::Options opt, Args_&&... args) : 
            UnknownExtractor<byrow_, selection_, true>(p, std::forward<Args_>(args)...), options(std::move(opt)) {}

        tatami::SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            auto output = this->parent->template run_sparse_extractor<byrow_>(i, vbuffer, ibuffer, options, this->work.get());

            // Need to adjust the indices.
            if (output.index) {
                if constexpr(selection_ == tatami::DimensionSelectionType::BLOCK) {
                    for (Index_ i = 0; i < output.number; ++i) {
                        ibuffer[i] = output.index[i] + this->block_start;
                    }
                    output.index = ibuffer;
                } else if constexpr(selection_ == tatami::DimensionSelectionType::INDEX) {
                    for (Index_ i = 0; i < output.number; ++i) {
                        ibuffer[i] = this->indices[output.index[i]];
                    }
                    output.index = ibuffer;
                }
            }

            return output;
        }

    private:
        tatami::Options options;
    };

public:
    std::unique_ptr<tatami::FullDenseExtractor<Value_, Index_> > dense_row(const tatami::Options& opt) const {
        return std::unique_ptr<tatami::FullDenseExtractor<Value_, Index_> >(new DenseUnknownExtractor<true, tatami::DimensionSelectionType::FULL>(this, opt));
    }

    std::unique_ptr<tatami::BlockDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const tatami::Options& opt) const {
        return std::unique_ptr<tatami::BlockDenseExtractor<Value_, Index_> >(new DenseUnknownExtractor<true, tatami::DimensionSelectionType::BLOCK>(this, opt, block_start, block_length));
    }

    std::unique_ptr<tatami::IndexDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const tatami::Options& opt) const {
        return std::unique_ptr<tatami::IndexDenseExtractor<Value_, Index_> >(new DenseUnknownExtractor<true, tatami::DimensionSelectionType::INDEX>(this, opt, std::move(indices)));
    }

    std::unique_ptr<tatami::FullDenseExtractor<Value_, Index_> > dense_column(const tatami::Options& opt) const {
        return std::unique_ptr<tatami::FullDenseExtractor<Value_, Index_> >(new DenseUnknownExtractor<false, tatami::DimensionSelectionType::FULL>(this, opt));
    }

    std::unique_ptr<tatami::BlockDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const tatami::Options& opt) const {
        return std::unique_ptr<tatami::BlockDenseExtractor<Value_, Index_> >(new DenseUnknownExtractor<false, tatami::DimensionSelectionType::BLOCK>(this, opt, block_start, block_length));
    }

    std::unique_ptr<tatami::IndexDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const tatami::Options& opt) const {
        return std::unique_ptr<tatami::IndexDenseExtractor<Value_, Index_> >(new DenseUnknownExtractor<false, tatami::DimensionSelectionType::INDEX>(this, opt, std::move(indices)));
    }

public:
    std::unique_ptr<tatami::FullSparseExtractor<Value_, Index_> > sparse_row(const tatami::Options& opt) const {
        return std::unique_ptr<tatami::FullSparseExtractor<Value_, Index_> >(new SparseUnknownExtractor<true, tatami::DimensionSelectionType::FULL>(this, opt));
    }

    std::unique_ptr<tatami::BlockSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_length, const tatami::Options& opt) const {
        return std::unique_ptr<tatami::BlockSparseExtractor<Value_, Index_> >(new SparseUnknownExtractor<true, tatami::DimensionSelectionType::BLOCK>(this, opt, block_start, block_length));
    }

    std::unique_ptr<tatami::IndexSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> indices, const tatami::Options& opt) const {
        return std::unique_ptr<tatami::IndexSparseExtractor<Value_, Index_> >(new SparseUnknownExtractor<true, tatami::DimensionSelectionType::INDEX>(this, opt, std::move(indices)));
    }

    std::unique_ptr<tatami::FullSparseExtractor<Value_, Index_> > sparse_column(const tatami::Options& opt) const {
        return std::unique_ptr<tatami::FullSparseExtractor<Value_, Index_> >(new SparseUnknownExtractor<false, tatami::DimensionSelectionType::FULL>(this, opt));
    }

    std::unique_ptr<tatami::BlockSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_length, const tatami::Options& opt) const {
        return std::unique_ptr<tatami::BlockSparseExtractor<Value_, Index_> >(new SparseUnknownExtractor<false, tatami::DimensionSelectionType::BLOCK>(this, opt, block_start, block_length));
    }

    std::unique_ptr<tatami::IndexSparseExtractor<Value_, Index_> > sparse_column(std::vector<Index_> indices, const tatami::Options& opt) const {
        return std::unique_ptr<tatami::IndexSparseExtractor<Value_, Index_> >(new SparseUnknownExtractor<false, tatami::DimensionSelectionType::INDEX>(this, opt, std::move(indices)));
    }
};

}

#endif
