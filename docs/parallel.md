# Enabling parallelization

By default, it is assumed that `tatami_r::UnknownMatrix` will only be used in a single-threaded context.
This is because the `UnknownMatrix` methods will call the R API, which is strictly single-threaded.
In fact, it's worse than that: some experimentation indicates that R code can only be executed on the main thread.
Otherwise, we get stack limit errors when R is called inside a worker, even with locking to enforce serial execution.
I suspect that some other R-managed process (an event loop, perhaps?) is always running on the main thread;
calling R from the worker will fail to block this process and lead to parallel execution between the worker and the main threads.

As a result, some work is required to allow `UnknownMatrix` objects to operate in a parallel context.
Here, we will be considering worker threads spawned via the `std::thread` library, where each worker may call methods of a `tatami_r::UnknownMatrix`.
We use the [**manticore**](https://github.com/tatami-inc/manticore) library to allow each worker thread to pass a function for execution on the main thread;
once execution is complete, the results are converted into something independent of R, and those are returned to the worker.
This ensures that the R runtime is safely evaluated without requiring all **Rcpp** code to be lifted out of the `tatami::Matrix`'s getter methods
(which would otherwise require a complete redesign of all applications that consume `tatami::Matrix` objects),

Enabling a **manticore** parallel context is simple.
Assume that we already have a `tatami::Matrix` object that _might_ contain a `UnknownMatrix`.
To enable safe parallel execution, we call `initialize()` before the parallel section and `listen()` afterwards.
This diverts the R code to the main thread for execution while allowing all other **tatami**-related code to run inside each worker.
Note that users must define the `TATAMI_R_PARALLELIZE_UNKNOWN` macro to expose the `executor()` function.

```cpp
#define TATAMI_R_PARALLELIZE_UNKNOWN

// Initializing the context with a manticore::Executor.
auto& mexec = tatami_r::executor();
mexec.initialize(num_threads);

// Spin up std::thread objects. You can actually use anything
// that respects <mutex>, <atomic> and <condition_variable>.
std::vector<std::thread> threads;
threads.reserve(num_threads);

for (int t = 0; t < num_threads; ++t) {
    threads.emplace_back([&]() -> void {
        // Do something with the tatami:Matrix here.
        auto wrk = my_matrix->dense_row();
        auto row_t = wrk->fetch(t);
        mexec.finish_thread();
    });
}

// Listen for requests from the worker threads.
mexec.listen();

// Join all threads once the work is done.
for (auto& th : threads) {
    th.join();
}
```

Even more simply, we provide the pre-packaged `tatami_r::parallelize()` method for use as a custom **tatami** parallelization scheme.
This requires a little bit of care to use, as the `TATAMI_CUSTOM_PARALLEL` macro has to be set before any includes of **tatami** (which is done automatically by including `tatami_r/tatami_r.hpp`).
We suggest using the following sequence of calls before any **tatami**-containing source file.

```cpp
// Must be defined before including tatami_r/parallelize.hpp
#define TATAMI_R_PARALLELIZE_UNKNOWN 
#include "tatami_r/parallelize.hpp"

// Defining the tatami parallelization scheme, which must 
// occur before including tatami or tatami_r itself.
#define TATAMI_CUSTOM_PARALLEL tatami_r::parallelize

#include "tatami_r/tatami_r.hpp"
```

Note that construction of the `UnknownMatrix` must always be performed on the main thread.
This is because construction of the unknown fallback involves some calls into the R runtime;
these are currently not protected from execution in worker contexts.
