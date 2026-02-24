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
#include "sanisizer/sanisizer.hpp"

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
 * @cond
 */
inline manticore::Executor* executor_ptr = NULL;
/**
 * @endcond
 */

/**
 * Retrieve a global `manticore::Executor` object for all **tatami_r** applications.
 * This function is only available if `TATAMI_R_PARALLELIZE_UNKNOWN` is defined.
 *
 * @return Reference to a global `manticore::Executor`.
 * If `set_executor()` was called with a non-`NULL` pointer, the provided instance will be used;
 * otherwise, a default instance will be instantiated.
 */
inline manticore::Executor& executor() {
    if (executor_ptr) {
        return *executor_ptr;
    } else {
        // In theory, this should end up resolving to a single instance, even across dynamically linked libraries:
        // https://stackoverflow.com/questions/52851239/local-static-variable-linkage-in-a-template-class-static-member-function
        // In practice, this doesn't seem to be the case on a Mac, requiring us to use `set_executor()`.
        static manticore::Executor mexec;
        return mexec;
    }
}

/**
 * Set a global `manticore::Executor` object for all **tatami_r** applications.
 * This function is only available if `TATAMI_R_PARALLELIZE_UNKNOWN` is defined.
 * Calling this function is occasionally necessary if `executor()` resolves to different instances of a `manticore::Executor` across different libraries.
 *
 * @param Pointer to a global `manticore::Executor`, or `NULL` to unset this pointer.
 */
inline void set_executor(manticore::Executor* ptr) {
    executor_ptr = ptr;
}

/**
 * @tparam Function_ Function to be executed.
 * @tparam Index_ Integer type for the task indices.
 *
 * @param fun Function to run in each thread.
 * This is a lambda that should accept three arguments:
 * - Integer containing the thread ID in `[0, threads)`. 
 * - Integer specifying the index of the first task to be executed in a thread.
 *   This will lie in `[0, tasks)`.
 * - Integer specifying the number of tasks to be executed in a thread.
 *   This will lie in `(0, tasks)`, i.e., it is always positive.
 * @param tasks Number of tasks to be executed.
 * This should be non-negative.
 * @param threads Number of threads to parallelize over.
 * This should be positive.
 *
 * This function is a drop-in replacement for `tatami::parallelize()`.
 * The series of integers from `[0, tasks)` is deterministically split into `K` non-overlapping non-empty contiguous ranges where `K <= threads`.
 * Each range is passed to `fun` for parallel execution with the standard `<thread>` library. 
 * Serialization can be achieved via `<mutex>` in most cases, or `manticore::Executor::run()` if the task must be performed on the main thread (see `executor()`).
 *
 * This function is only available if `TATAMI_R_PARALLELIZE_UNKNOWN` is defined.
 *
 * @return The number of workers (`K`) that were actually used.
 * `K` is guaranteed to be no greater than `threads` (or 1, if the latter is not positive).
 * `fun()` will have been called once for each of the thread IDs `[0, ..., K - 1]`.
 */ 
template<class Function_, class Index_>
int parallelize(const Function_ fun, const Index_ tasks, int threads) {
    if (tasks == 0) {
        return 0;
    }

    if (threads <= 1 || tasks == 1) {
        fun(0, 0, tasks);
        return 1;
    }

    Index_ tasks_per_worker = tasks / threads;
    int remainder = tasks % threads;
    if (tasks_per_worker == 0) {
        tasks_per_worker = 1; 
        remainder = 0;
        threads = tasks;
    }

    auto& mexec = executor();
    mexec.initialize(threads, "failed to execute R command");

    std::vector<std::thread> runners;
    sanisizer::reserve(runners, threads);
    auto errors = sanisizer::create<std::vector<std::exception_ptr> >(threads);

    Index_ start = 0;
    for (int w = 0; w < threads; ++w) {
        Index_ length = tasks_per_worker + (w < remainder);

        runners.emplace_back(
            [&](const int id, const Index_ s, const Index_ l) -> void {
                try {
                    fun(id, s, l);
                } catch (...) {
                    errors[id] = std::current_exception();
                }
                mexec.finish_thread();
            },
            w,
            start,
            length
        );

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

    return threads;
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
