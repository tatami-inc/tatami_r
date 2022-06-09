#ifndef RATICATE_UNKNOWNMATRIX_HPP
#define RATICATE_UNKNOWNMATRIX_HPP

#include "Rcpp.h"
#include "tatami/tatami.hpp"
#include "SimpleMatrix.hpp"
#include "SparseArraySeed.hpp"
#include <vector>
#include <memory>
#include <string>
#include <stdexcept>

namespace raticate {

template<typename Data, typename Index>
class UnknownMatrix : public tatami::Matrix<Data, Index> {
public: 
    UnknownMatrix(Rcpp::RObject seed) :
        original_seed(seed),
        delayed_env(Rcpp::Environment::namespace_env("DelayedArray")),
        dense_extractor(delayed_env["extract_array"]),
        sparse_extractor(delayed_env["extract_sparse_array"])
    {
        // We assume the constructor is already wrapped in a serial section by
        // the caller, so we won't bother adding the various omp critical
        // statements here. I'm also not sure that the operations in the
        // initialization list are thread-safe.

        {
            auto base = Rcpp::Environment::base_env();
            Rcpp::Function fun = base["dim"];
            Rcpp::RObject output = fun(seed);
            if (output.sexp_type() != INTSXP) {
                throw std::runtime_error("'dims' should return an integer vector");
            }
            Rcpp::IntegerVector dims(output);
            if (dims.size() != 2 || dims[0] < 0 || dims[1] < 0) {
                throw std::runtime_error("'dims' should contain two non-negative integers");
            }
            nrow_ = dims[0];
            ncol_ = dims[1];
        }

        {
            Rcpp::Function fun = delayed_env["is_sparse"];
            Rcpp::LogicalVector sparse = fun(seed);
            if (sparse.size() != 1) {
                throw std::runtime_error("'is_sparse' should return a logical vector of length 1");
            }
            sparse_ = (sparse[0] != 0);
        }

        {
            Rcpp::Function fun = delayed_env["chunkdim"];
            Rcpp::RObject output = fun(seed);
            needs_chunks = !output.isNULL();
            if (needs_chunks) {
                Rcpp::IntegerVector chunks(output);
                if (chunks.size() != 2 || chunks[0] < 0 || chunks[1] < 0) {
                    throw std::runtime_error("'chunks' should contain two non-negative integers");
                }
                chunk_nrow = chunks[0];
                chunk_ncol = chunks[1];
            }
        }

        {
            Rcpp::Function fun = delayed_env["colAutoGrid"];
            Rcpp::RObject output = fun(seed);
            Rcpp::IntegerVector spacing = output.slot("spacings");
            if (spacing.size() != 2 || spacing[1] < 0) {
                throw std::runtime_error("'spacings' slot of 'colAutoGrid' output should contain two non-negative integers");
            }
            block_ncol = spacing[1];
        }

        {
            Rcpp::Function fun = delayed_env["rowAutoGrid"];
            Rcpp::RObject output = fun(seed);
            Rcpp::IntegerVector spacing = output.slot("spacings");
            if (spacing.size() != 2 || spacing[0] < 0) {
                throw std::runtime_error("'spacings' slot of 'rowAutoGrid' output should contain two non-negative integers");
            }
            block_nrow = spacing[0];
        }
    }

public:
    size_t nrow() const {
        return nrow_;
    }

    size_t ncol() const {
        return ncol_;
    }

    bool sparse() const {
        return sparse_;
    }

    bool prefer_rows() const {
        // All of the individual extract_array outputs are effectively column-major.
        return false;
    }

public:
    struct UnknownWorkspace : public tatami::Workspace {
        UnknownWorkspace(bool r = true) : byrow(r) {}
        bool byrow;

        size_t primary_block_start, primary_block_end;
        size_t secondary_chunk_start, secondary_chunk_end;

        std::shared_ptr<tatami::Matrix<Data, Index> > buffer = nullptr;
        std::shared_ptr<tatami::Workspace> bufwork = nullptr;

        Rcpp::RObject contents;
    };

    std::shared_ptr<tatami::Workspace> new_workspace(bool row) const { 
        std::shared_ptr<tatami::Workspace> output;

#ifndef RATICATE_RCPP_PARALLEL_LOCK        
        #pragma omp critical(RATICATE_RCPP_CRITICAL_NAME)
        {
#else
        RATICATE_RCPP_PARALLEL_LOCK([&]() -> void {
#endif

            std::cout << "Creating..." << std::endl;            
            output.reset(new UnknownWorkspace(row));
            std::cout << "READY!" << std::endl;            
                     
#ifndef RATICATE_RCPP_PARALLEL_LOCK        
        }
#else
        });
#endif

        return output;
    }

private:
    static Rcpp::RObject create_index_vector(size_t first, size_t last, size_t max) const {
        if (first != 0 || last != max) {
            Rcpp::IntegerVector alt(last - first);
            std::iota(alt.begin(), alt.end(), first + 1); // 1-based.
            return alt;
        } else {
            return R_NilValue;
        }
    }

    template<bool byrow>
    Rcpp::List create_quick_indices(size_t i, size_t first, size_t last) const {
        Rcpp::List indices(2);
        indices[(byrow ? 0 : 1)] = Rcpp::IntegerVector::create(i + 1);
        indices[(byrow ? 1 : 0)] = create_index_vector(first, last, (byrow ? ncol_ : nrow_));
        return indices;
    }

    std::pair<size_t, size_t> round_indices(size_t first, size_t last, size_t interval, size_t max) const {
        size_t new_first = (first / interval) * interval;
        size_t new_last = std::min(
            max, 
            (last ? 
                ((last - 1) / interval + 1) * interval // i.e., ceil(last/interval) * interval.
                : 0 
            )
        );
        return std::make_pair(new_first, new_last);
    }

    template<bool byrow>
    Rcpp::List create_rounded_indices(size_t i, size_t first, size_t last, UnknownWorkspace* work) const {
        Rcpp::List indices(2);
        if constexpr(byrow) {
            auto row_rounded = round_indices(i, i + 1, block_nrow, nrow_);
            indices[0] = create_index_vector(row_rounded.first, row_rounded.second, nrow_);
            work->primary_block_start = row_rounded.first;
            work->primary_block_end = row_rounded.second;

            auto col_rounded = (needs_chunks ? round_indices(first, last, chunk_ncol, ncol_) : std::make_pair(first, last));
            indices[1] = create_index_vector(col_rounded.first, col_rounded.second, ncol_);
            work->secondary_chunk_start = col_rounded.first;
            work->secondary_chunk_end = col_rounded.second;

        } else {
            auto row_rounded = (needs_chunks ? round_indices(first, last, chunk_nrow, nrow_) : std::make_pair(first, last));
            indices[0] = create_index_vector(row_rounded.first, row_rounded.second, nrow_);
            work->secondary_chunk_start = row_rounded.first;
            work->secondary_chunk_end = row_rounded.second;

            auto col_rounded = round_indices(i, i + 1, block_ncol, ncol_);
            indices[1] = create_index_vector(col_rounded.first, col_rounded.second, ncol_);
            work->primary_block_start = col_rounded.first;
            work->primary_block_end = col_rounded.second;
        }
        return indices;
    }

    bool needs_reset(size_t i, size_t first, size_t last, const UnknownWorkspace* work) const {
        bool reset = true;
        if (work->buffer != nullptr) {
            if (i >= work->primary_block_start && i < work->primary_block_end) {
                if (first >= work->secondary_chunk_start && last <= work->secondary_chunk_end) {
                    reset = false;
                }
            }
        }
        return reset;
    }

public:
    const Data* row(size_t r, Data* buffer, size_t first, size_t last, tatami::Workspace* work=nullptr) const {
        if (work == NULL) {
        } else {
            buffered_dense_extractor<true>(r, buffer, first, last, work);
        }
        return buffer;
    }

    const Data* column(size_t c, Data* buffer, size_t first, size_t last, tatami::Workspace* work=nullptr) const {
        if (work == NULL) {
            quick_dense_extractor<false>(c, buffer, first, last);
        } else {
            buffered_dense_extractor<false>(c, buffer, first, last, work);
        }
        return buffer;
    } 

public:
    /**
     * @cond
     */




    /**
     * @endcond
     */



    struct MainWorker {

    };
    

private:
    template<bool byrow>
    static void quick_dense_extractor_raw(size_t i, Data* buffer, size_t first, size_t last, const Rcpp::RObject* original_ptr, const Rcpp::Function* dense_extractor) {
        auto indices = create_quick_indices<byrow>(i, first, last);
        Rcpp::RObject val0 = (*dense_extractor)(*original_seed, indices);
        if (val0.sexp_type() == LGLSXP) {
            Rcpp::LogicalVector val(val0);
            std::copy(val.begin(), val.end(), buffer);
        } else if (val0.sexp_type() == INTSXP) {
            Rcpp::IntegerVector val(val0);
            std::copy(val.begin(), val.end(), buffer);
        } else {
            Rcpp::NumericVector val(val0);
            std::copy(val.begin(), val.end(), buffer);
        }
    }

    template<bool byrow>
    void quick_dense_extractor(size_t i, Data* buffer, size_t first, size_t last) const {
#ifndef RATICATE_RCPP_PARALLEL_LOCK
        quick_dense_extractor_raw<byrow>(i, buffer, first, last, &original_seed, &dense_extractor);
#else
        RATICATE_RCPP_PARALLEL([&]() -> void {
            executor().set<byrow>(i, buffer, first, last, &original_seed, &dense_extractor);
        });
#endif
    }

    template<bool byrow>
    static void buffered_dense_extractor_raw(size_t i, Data* buffer, size_t first, size_t last, UnknownWorkspace* work, const Rcpp::RObject* original_ptr, const Rcpp::Function* dense_extractor) {
        auto indices = create_rounded_indices<byrow>(i, first, last, work);
        Rcpp::RObject val0 = (*dense_extractor)(original_seed, indices);
        auto parsed = parse_simple_matrix<Data, Index>(val0);
        work->buffer = parsed.matrix;
        work->contents = parsed.contents;
        work->bufwork = (work->buffer)->new_workspace(byrow);
    }

    template<bool byrow>
    void buffered_dense_extractor(size_t i, Data* buffer, size_t first, size_t last, tatami::Workspace* work0) const {
        UnknownWorkspace* work = static_cast<UnknownWorkspace*>(work0);
        if (work->byrow != byrow) {
            throw std::runtime_error("workspace should have been generated with 'row=" + std::to_string(byrow) + "'");
        }

        if (needs_reset(i, first, last, work)) {
#ifndef RATICATE_RCPP_PARALLEL_LOCK
            buffered_dense_extractor_raw<byrow>(i, buffer, first, last, work, &original_seed, &dense_extractor);
#else
            RATICATE_RCPP_PARALLEL([&]() -> void {
                executor().set<byrow>(i, buffer, first, last, work, &original_seed, &dense_extractor);
            });
#endif
        }

        i -= work->primary_block_start;
        first -= work->secondary_chunk_start;
        last -= work->secondary_chunk_start;
        if constexpr(byrow) {
            (work->buffer)->row_copy(i, buffer, first, last, (work->bufwork).get());
        } else {
            (work->buffer)->column_copy(i, buffer, first, last, (work->bufwork).get());
        }
    }

public:
    virtual tatami::SparseRange<Data, Index> sparse_row(size_t r, Data* vbuffer, Index* ibuffer, size_t first, size_t last, tatami::Workspace* work=nullptr, bool sorted=true) const {
        if (sparse_) {
            if (work == NULL) {
                return quick_sparse_extractor<true>(r, vbuffer, ibuffer, first, last, sorted);
            } else {
                return buffered_sparse_extractor<true>(r, vbuffer, ibuffer, first, last, work, sorted);
            }
        } else {
            return tatami::Matrix<Data, Index>::sparse_row(r, vbuffer, ibuffer, first, last, work, sorted);
        }
    }

    virtual tatami::SparseRange<Data, Index> sparse_column(size_t c, Data* vbuffer, Index* ibuffer, size_t first, size_t last, tatami::Workspace* work=nullptr, bool sorted=true) const {
        if (sparse_) {
            if (work == NULL) {
                return quick_sparse_extractor<false>(c, vbuffer, ibuffer, first, last, sorted);
            } else {
                return buffered_sparse_extractor<false>(c, vbuffer, ibuffer, first, last, work, sorted);
            }
        } else {
            return tatami::Matrix<Data, Index>::sparse_column(c, vbuffer, ibuffer, first, last, work, sorted);
        }
    }

private:
    template<bool byrow>
    tatami::SparseRange<Data, Index> quick_sparse_extractor(size_t i, Data* vbuffer, Index* ibuffer, size_t first, size_t last, bool sorted) const {
        size_t n = 0;

#ifndef RATICATE_RCPP_PARALLEL_LOCK
        #pragma omp critical(RATICATE_RCPP_CRITICAL_NAME)
        {
#else
        RATICATE_RCPP_PARALLEL_LOCK([&]() -> void {
#endif

            auto indices = create_quick_indices<byrow>(i, first, last);
            Rcpp::RObject val0 = sparse_extractor(original_seed, indices);

            {
                Rcpp::IntegerMatrix indices = val0.slot("nzindex");
                n = indices.rows();
                auto idx = indices.column(byrow ? 1 : 0);
                auto icopy = ibuffer;
                for (auto ix : idx) {
                    *icopy = ix + first - 1; // 0-based indices.
                    ++icopy;
                }
            }

            Rcpp::RObject data = val0.slot("nzdata");
            if (data.sexp_type() == LGLSXP) {
                Rcpp::LogicalVector val(data);
                std::copy(val.begin(), val.end(), vbuffer);
            } else if (data.sexp_type() == INTSXP) {
                Rcpp::IntegerVector val(data);
                std::copy(val.begin(), val.end(), vbuffer);
            } else {
                Rcpp::NumericVector val(data);
                std::copy(val.begin(), val.end(), vbuffer);
            }

#ifndef RATICATE_RCPP_PARALLEL_LOCK        
        }
#else
        });
#endif

        if (sorted && !std::is_sorted(ibuffer, ibuffer + n)) {
            // TODO: use an in-place sort?
            std::vector<std::pair<Data, Index> > holding;
            holding.reserve(n);
            for (size_t ix = 0; ix < n; ++ix) {
                holding.emplace_back(vbuffer[ix], ibuffer[ix]);
            }
            std::sort(holding.begin(), holding.end());
            for (size_t ix = 0; ix < n; ++ix) {
                vbuffer[ix] = holding[ix].first;
                ibuffer[ix] = holding[ix].second;
            }
        }

        return tatami::SparseRange<Data, Index>(n, vbuffer, ibuffer);
    }

    template<bool byrow>
    tatami::SparseRange<Data, Index> buffered_sparse_extractor(size_t i, Data* vbuffer, Index* ibuffer, size_t first, size_t last, tatami::Workspace* work0, bool sorted) const {
        UnknownWorkspace* work = static_cast<UnknownWorkspace*>(work0);
        if (work->byrow != byrow) {
            throw std::runtime_error("workspace should have been generated with 'row=" + std::to_string(byrow) + "'");
        }

        if (needs_reset(i, first, last, work)) {
#ifndef RATICATE_RCPP_PARALLEL_LOCK
            #pragma omp critical(RATICATE_RCPP_CRITICAL_NAME)
            {
#else
            RATICATE_RCPP_PARALLEL_LOCK([&]() -> void {
#endif

                auto indices = create_rounded_indices<byrow>(i, first, last, work);
                auto val0 = sparse_extractor(original_seed, indices);
                auto parsed = parse_SparseArraySeed<Data, Index>(val0);
                work->buffer = parsed.matrix;
                work->contents = parsed.contents;
                work->bufwork = (work->buffer)->new_workspace(byrow);

#ifndef RATICATE_RCPP_PARALLEL_LOCK        
            }
#else
            });
#endif
        }

        i -= work->primary_block_start;
        first -= work->secondary_chunk_start;
        last -= work->secondary_chunk_start;

        tatami::SparseRange<Data, Index> output;
        if constexpr(byrow) {
            output = (work->buffer)->sparse_row_copy(i, vbuffer, ibuffer, first, last, tatami::SPARSE_COPY_BOTH, (work->bufwork).get(), sorted);
        } else {
            output = (work->buffer)->sparse_column_copy(i, vbuffer, ibuffer, first, last, tatami::SPARSE_COPY_BOTH, (work->bufwork).get(), sorted);
        }

        // Need to adjust the indices.
        for (size_t i = 0; i < output.number; ++i) {
            ibuffer[i] += work->secondary_chunk_start;
        }

        return output;
    }

private:
    size_t nrow_, ncol_;
    bool sparse_;

    bool needs_chunks;
    size_t chunk_nrow, chunk_ncol;

    size_t block_nrow, block_ncol;

public:
    friend struct Executor {
        bool sparse;
        bool buffered;
        bool byrow;

        size_t index, first, last;
        Data* dbuffer;
        Index* ibuffer;
        UnknownMatrix::Workspace* work;

        const Rcpp::RObject* original_ptr;
        const Rcpp::Function* extractor;

        bool ready = false;
        bool finished = false;

        template<bool B>
        void set(size_t i, size_t f, size_t l, const Rcpp::RObject* original_seed) {
            byrow = B;
            index = i;
            first = f;
            last = l;
            original_ptr = original_seed;
            ready = true;
            finished = false;
        }

        template<bool B>
        void set(size_t i, Data* buffer, size_t f, size_t l, const Rcpp::RObject* original_seed, const Rcpp:Function* dense_extractor) {
            set(i, f, l, original_seed);
            sparse = false;
            buffered = false;
            dbuffer = buffer;
            extractor = dense_extractor;
        }

        template<bool B>
        void set(size_t i, Data* buffer, size_t f, size_t l, UnknownMatrix::Workspace* w, const Rcpp::RObject* original_seed, const Rcpp:Function* dense_extractor) {
            set(i, buffer, f, l, original_seed, dense_extractor);
            buffered = true;
            work = w;
        }

        void harvest() {
            if (!sparse) {
                if (buffered) {
                    if (byrow) {
                        buffered_dense_extractor_raw<true>(i, first, last, work, original_ptr, dense_extractor);
                    } else {
                        buffered_dense_extractor_raw<false>(i, first, last, work, original_ptr, dense_extractor);
                    }
                } else {
                    if (byrow) {
                        quick_dense_extractor_raw<true>(i, first, last, original_ptr, dense_extractor);
                    } else {
                        quick_dense_extractor_raw<false>(i, first, last, original_ptr, dense_extractor);
                    }
                }
            }
            finished = true;
        }
    };

    static Executor& executor() {
        static Executor ex;
        return ex;
    }

private:
    Rcpp::RObject original_seed;
    Rcpp::Environment delayed_env;
    Rcpp::Function dense_extractor, sparse_extractor;
};

}

#endif
