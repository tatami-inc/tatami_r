#ifndef RATICATE_PARALLEL_HPP
#define RATICATE_PARALLEL_HPP

namespace raticate {

class Parallelizer {
    static std::mutex rcpp_lock;
    static std::condition_variable cv;

    static bool ready;
    static bool processed;

    static size_t i;
    static size_t first;
    static size_t last;
    static const Rcpp::RObject* seed_ptr;

    static bool quick;
    static bool sparse;

    // Dense case.
    static void* buffer;

    // Sparse case.
    static void* ibuffer;
    static void* vbuffer;
    static bool sorted;
};

template<class Function, typename Data, typename Index>
void parallelize(Function fun, size_t njobs, size_t nthreads) {
    size_t jobs_per_worker = std::ceil(static_cast<double>(njobs) / nthreads);
    size_t start = 0;
    std::vector<std::thread> jobs;
    
    for (size_t w = 0; w < nthreads; ++w) {
        size_t end = std::min(n, start + jobs_per_worker);
        if (start >= end) {
            break;
        }
        jobs.emplace_back(f, start, end);
        start += jobs_per_worker;
    }

    // Main thread waits for all jobs to exhaust themselves.
    while (1) {
        std::unique_lock lk(Parallelizer::rcpp_lock);

        Parallelizer::cv.wait(lk, []{ return Parallelizer::ready; });
        Parallelizer::ready = false;
 
        if (Parallelizer::quick) {
            Rcpp::List indices;
            if (Parallelizer::byrow) {
                indices = create_quick_indices<true>(Parallelizer::i, Parallelizer::first, Parallelizer::last);
            } else {
                indices = create_quick_indices<false>(Parallelizer::i, Parallelizer::first, Parallelizer::last);
            }

            if (Parallelizer::dense) {
                Rcpp::RObject val0 = dense_extractor(*Parallelizer::seed_ptr, indices);
                Data* buffer = static_cast<Data*>(Parallelizer::buffer);
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
        }
    
        Parallelizer::processed = true;

        // Manual unlocking is done before notifying, to avoid waking up
        // the waiting thread only to block again (see notify_one for details)
        lk.unlock();
        cv.notify_one();
    }

    for (auto& job : jobs) {
        job.join();
    }

}

template<class Function, class Data>
void quick_dense_locks(size_t i, size_t first, size_t last, Data* buffer, const Rcpp::RObject& obj) {
    {
        std::lock_guard lk(Parallelizer::rcpp_lock);
        Parallelizer::ready = true;
        Paralellizer::processed = false;
        Parallelizer::i = i;
        Parallelizer::first = first;
        Parallelizer::last = last;
        Parallelizer::buffer = buffer;
        Parallelizer::seed_ptr = &obj;
        Parallelizer::quick = true;
        Parallelizer::sparse = false;
    }

    Parallelizer::cv.notify_one();

    {
        std::unique_lock<std::mutex> lk(Parallelizer::rcpp_lock);
        Parallelizer::cv.wait(lk, []{ return Parallelizer::processed; });
    }
}

}

#endif
