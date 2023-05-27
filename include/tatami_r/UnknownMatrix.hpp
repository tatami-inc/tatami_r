#ifndef TATAMI_R_UNKNOWNMATRIX_HPP
#define TATAMI_R_UNKNOWNMATRIX_HPP

#include "Rcpp.h"
#include "tatami/tatami.hpp"
#include "SimpleMatrix.hpp"
#include "SparseArraySeed.hpp"
#include <vector>
#include <memory>
#include <string>
#include <stdexcept>
#include <unordered_map>

#include "parallelize.hpp"

namespace tatami_r {

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
        dense_extractor(delayed_env["extract_array"]),
        sparse_extractor(delayed_env["OLD_extract_sparse_array"])
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
            Rcpp::Function fun = delayed_env["chunkdim"];
            Rcpp::RObject output = fun(seed);
            if (output == R_NilValue) {
                chunk_nrow = 1;
                chunk_ncol = 1;
            } else {
                Rcpp::IntegerVector chunks(output);
                if (chunks.size() != 2 || chunks[0] < 0 || chunks[1] < 0) {
                    auto ctype = get_class_name(original_seed);
                    throw std::runtime_error("'chunkdim(<" + ctype + ">)' should return two non-negative integers");
                }
                chunk_nrow = chunks[0];
                chunk_ncol = chunks[1];
            }
        }

        cache_size = cache;
        if (cache_size == -1) {
            Rcpp::Function fun = delayed_env["getAutoBlockSize"];
            Rcpp::NumericVector output = fun();
            if (output.size() != 1 || output[0] < 0) {
                throw std::runtime_error("'getAutoBlockSize()' should return a non-negative number of bytes");
            }
            cache_size = output[0];
        }

        auto chunks_per_row = static_cast<double>(internal_ncol) / chunk_ncol;
        auto chunks_per_col = static_cast<double>(internal_nrow) / chunk_nrow;
        internal_prefer_rows = chunks_per_row <= chunks_per_col;
    }

private:
    Index_ internal_nrow, internal_ncol;
    bool internal_sparse, internal_prefer_row;

    size_t cache_size;
    Index_ chunk_nrow, chunk_ncol;

    Rcpp::RObject original_seed;
    Rcpp::Environment delayed_env;
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

    bool prefer_rows() const {
        return internal_prefer_rows;
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

        Index_ block_size;

        // No oracle, we just go for blocks.
        Index_ primary_block_start, primary_block_len;

        // With an oracle, we try to make some better decisions.
        tatami::OracleStream<Index_> prediction_stream;
        std::unordered_map<Index_, Index_> predictors;
        std::vector<Index_> unique_predictions;
        size_t predictions_made = 0, predictions_used = 0;
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
            auto row_rounded = round_indices(i, work->block_size, internal_nrow);
            work->primary_block_start = row_rounded.first;
            work->primary_block_len = row_rounded.second;
            indices[0] = create_consecutive_indices(row_rounded.first, row_rounded.second);
            indices[1] = work->secondary_indices;
        } else {
            auto col_rounded = round_indices(i, work->block_size, internal_ncol);
            work->primary_block_start = col_rounded.first;
            work->primary_block_len = col_rounded.second;
            indices[0] = work->secondary_indices;
            indices[1] = create_consecutive_indices(col_rounded.first, col_rounded.second);
        }
        return indices;
    }

    template<bool byrow_, bool sparse_>
    Rcpp::List create_next_indices(Index_ i, Workspace<sparse_>* work) const {
        work->unique_predictions.clear();
        work->predictors.clear();
        bool recycle = true;

        work->predictions_used = 0;
        size_t max_predictions = work->block_size * 10;
        auto& i = work->predictions_made;
        for (i = 0; i < max_predictions && work->unique_predictions.size() < work->block_size; ++i) {
            Index_ current;
            if (!work->prediction_stream.next(current)) {
                break;
            }

            auto it = work->predictors.find(current);
            if (it == work->predictors.end()) {
                work->predictors[current] = 0;
                unique_predictions.push_back(current);
            }
        }

        if (!std::is_sorted(unique_predictions.begin(), unique_predictions.end())) {
            std::sort(unique_predictions.begin(), unique_predictions.end());
        }
        Index_ counter = 0;
        for (auto x : unique_predictions) {
            work->predictors[x] = counter;
            ++counter;
        }
        Rcpp::IntegerVector primary_indices(unique_predictions.begin(), unique_predictions.end());
        for (auto& x : primary_indices) { // get to 1-based indices.
            ++x;
        }

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
        size_t parsed_primary = (byrow_ ? parsed->nrow() : parsed->ncol());
        size_t parsed_secondary = (byrow_ ? parsed->ncol() : parsed->nrow());

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
                run_dense_extractor<byrow_>(i, options, work, [](Index_ i2, Workspace<false>* work2) -> Rcpp::List { 
                    return create_next_indices<byrow_>(i2, work2);
                });
            }
            i = work->predictors.find(i)->second;
            ++work->predictions_used;

        } else {
            if (needs_reset(i, work)) {
                run_dense_extractor<byrow_>(i, options, work, [](Index_ i2, Workspace<false>* work2) -> Rcpp::List { 
                    return create_rounded_indices<byrow_>(i2, work2);
                });
            }
            i -= work->primary_block_start;
        }

        // Forcing a copy to avoid a possible reference to a transient buffer inside 'work->buffer'.
        return work->bufextractor->fetch_copy(i, buffer);
    }

    template<bool byrow_, class Function_>
    void run_sparse_extractor(Index_ i, Value_* buffer, const tatami::Options& options, Workspace<false>* work, Function_ index_creator) const {
#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
        // This involves some Rcpp initializations, so we lock it just in case.
        auto& mexec = executor();
        mexec.run([&]() -> void {
#endif

        auto indices = index_creator(i, work);

        if (internal_sparse) {
            auto val0 = sparse_extractor(original_seed, indices);
            auto parsed = parse_SparseArraySeed<Value_, Index_>(val0);
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

        check_buffered_dims<byrow_, false>(work->buffer.get(), work);
        work->bufextractor = tatami::new_extractor<byrow_, false>(work->buffer.get(), options);
    }

    template<bool byrow_>
    SparseRange<Value_, Index_> run_sparse_extractor(Index_ i, Value_* vbuffer, Index_* ibuffer, const tatami::Options& options, Workspace<true>* work) const {
        if (work->prediction_stream.active()) {
            if (work->predictions_used == work->predictions_made) {
                run_sparse_extractor<byrow_>(i, options, work, [](Index_ i2, Workspace<false>* work2) -> Rcpp::List { 
                    return create_next_indices<byrow_>(i2, work2);
                });
            }
            i = work->predictors.find(i)->second;
            ++work->predictions_used;

        } else {
            if (needs_reset(i, work)) {
                run_sparse_extractor<byrow_>(i, options, work, [](Index_ i2, Workspace<false>* work2) -> Rcpp::List { 
                    return create_rounded_indices<byrow_>(i2, work2);
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

        static void define_blocks(const UnknownMatrix<Value_, Index_>* parent, Index_ len, Workspace<sparse_>& work) {
            double cache_elements = static_cast<double>(parent->cache_size) / (static_cast<double>(len) * static_cast<double>(sizeof(Value_)));
            double chunk_dim = byrow_ ? parent->chunk_nrow : parent->chunk_ncol;
            work->block_size = std::max(1.0, std::floor(cache_elements / chunk_dim)) * chunk_dim;
        }

        UnknownExtractor(const UnknownMatrix<Value_, Index_>* p) : parent(p) { 
            if constexpr(selection_ == tatami::DimensionSelectionType::FULL) {
                this->full_length = byrow_ ? p->internal_ncol : p->internal_nrow;
                work.reset(setup_workspace(this->full_length));
            }
        }

        UnknownExtractor(const UnknownMatrix<Value_, Index_>* p, Index_ start, Index_ length) : parent(p) {
            if constexpr(selection_ == tatami::DimensionSelectionType::BLOCK) {
                this->block_start = start;
                this->block_length = length;
                work.reset(setup_workspace(start, length));
            }
        }

        UnknownExtractor(const UnknownMatrix<Value_, Index_>* p, std::vector<Index_> idx) : parent(p) {
            if constexpr(selection_ == tatami::DimensionSelectionType::INDEX) {
                indices = std::move(idx);
                this->index_length = indices.size();
                work.reset(setup_workspace(indices));
            }
        }

        const Index_* index_start() const {
            if constexpr(selection_ == tatami::DimensionSelectionType::INDEX) {
                return indices.data();
            } else {
                return NULL;
            }
        }

        void set_oracle(std::unique_ptr<tatami::Oracle<Index_> >) {
            // No-op.
            return;
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
    template<bool byrow, tatami::DimensionSelectionType selection_>
    struct SparseUnknownExtractor : public UnknownExtractor<byrow, selection_, true> {
        template<typename ... Args_>
        SparseUnknownExtractor(const UnknownMatrix<Value_, Index_>* p, tatami::Options opt, Args_&&... args) : 
            UnknownExtractor<byrow, selection_, true>(p, std::forward<Args_>(args)...), options(std::move(opt)) {}

        tatami::SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            return this->parent->template run_sparse_extractor<byrow_>(i, vbuffer, ibuffer, options, this->work.get());
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
