#ifndef TATAMI_R_PARALLELIZE_HPP
#define TATAMI_R_PARALLELIZE_HPP

/**
 * @cond
 */
#ifdef TATAMI_R_PARALLELIZE_UNKNOWN
/**
 * @endcond
 */

#include "manticore/manticore.hpp"
#include <thread>
#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>

/**
 * @file parallelize.hpp
 *
 * @brief Safely parallelize for unknown matrices.
 */

namespace tatami_r {

/**
 * Retrieve a global `manticore::Executor` object for all **tatami_r** uses.
 * This function is only available if `TATAMI_R_PARALLELIZE_UNKNOWN` is defined.
 *
 * @return Reference to a global `manticore::Executor`.
 */
inline manticore::Executor& executor() {
    // This should end up resolving to a single instance, even across dynamically linked libraries:
    // https://stackoverflow.com/questions/52851239/local-static-variable-linkage-in-a-template-class-static-member-function
    static manticore::Executor mexec;
    return mexec;
}

/**
 * @tparam Function_ Function to be executed.
 * @tparam Index_ Integer type for the job indices.
 *
 * @param fun Function to run in each thread.
 * This is a lambda that should accept three arguments:
 * - Integer containing the thread ID.
 * - Integer specifying the index of the first job to be executed in a thread.
 * - Integer specifying the number of jobs to be executed in a thread.
 * @param njobs Number of jobs to be executed.
 * @param nthreads Number of threads to parallelize over.
 *
 * This function is a drop-in replacement for `tatami::parallelize()`.
 * The series of integers from 0 to `njobs - 1` is split into `nthreads` contiguous ranges.
 * Each range is used as input to `fun` within the corresponding thread.
 * It is assumed that the execution of any given job is independent of the next.
 *
 * This function is only available if `TATAMI_R_PARALLELIZE_UNKNOWN` is defined.
 */ 
template<class Function_, class Index_>
void parallelize(Function_ fun, Index_ njobs, int nthreads) {
    if (njobs == 0) {
        return;
    }

    if (nthreads <= 1 || njobs == 1) {
        fun(0, 0, njobs);
        return;
    }

    Index_ jobs_per_worker = njobs / nthreads;
    int remainder = njobs % nthreads;
    if (jobs_per_worker == 0) {
        jobs_per_worker = 1; 
        remainder = 0;
        nthreads = njobs;
    }

    auto& mexec = executor();
    mexec.initialize(nthreads, "failed to execute R command");

    std::vector<std::thread> runners;
    runners.reserve(nthreads);
    std::vector<std::exception_ptr> errors(nthreads);

    Index_ start = 0;
    for (int w = 0; w < nthreads; ++w) {
        Index_ length = jobs_per_worker + (w < remainder);

        runners.emplace_back([&](int id, Index_ s, Index_ l) {
            try {
                fun(id, s, l);
            } catch (...) {
                errors[id] = std::current_exception();
            }
            mexec.finish_thread();
        }, w, start, length);

        start += length;
    }

    mexec.listen();
    for (auto& x : runners) {
        x.join();
    }

    for (const auto& err : errors) {
        if (err) {
            std::rethrow_exception(err);
        }
    }
}

}

/**
 * @cond
 */
#endif
/**
 * @endcond
 */

#endif
