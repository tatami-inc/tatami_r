# Enabling parallelization

## Overview

By default, we assume that the `tatami_r::UnknownMatrix` will only be used in a single-threaded context.
This is because the `UnknownMatrix` methods will call the R API, which is strictly single-threaded.
In fact, it's worse than that: some experimentation indicates that R code can only be executed on the main thread.
Otherwise, we get stack limit errors when R is called inside a worker, even with locking to enforce serial execution.
I suspect that some other R-managed process (an event loop, perhaps?) is always running on the main thread;
calling R from the worker will fail to block this process and lead to parallel execution between the worker and the main threads.

That said, it is possible to use `UnknownMatrix` objects in a parallel context with some effort.
Using the [**manticore**](https://github.com/tatami-inc/manticore) library, each worker thread can pass a function for execution on the main thread.
Specifically, multiple workers can request execution of the various `extract_*_array()` functions on the main thread;
once this is done, the extracted data is converted into something independent of R and returned to the worker.
Our assumption is that most of each worker's time is spent computing on the extracted data rather than waiting on the main thread,
such that parallelization still offers some performance improvement.

## Parallelizing matrix iterations

We provide a `tatami_r::parallelize()` function that handles all of the **manticore**-related boilerplate.
The `UnknownMatrix` methods will safely run within the lambda passed to this function.

```cpp
// Must be defined before including tatami_r/parallelize.hpp
#define TATAMI_R_PARALLELIZE_UNKNOWN 
#include "tatami_r/parallelize.hpp"

tatami_r::parallelize([&](size_t thread_id, int start, int len) -> void {
    // Do something with the UnknownMatrix.
    auto ext = ptr->dense_row();
    std::vector<double> buffer(ptr->ncol());
    for (int r = start, end = start + len; start < end; ++r) {
        auto out = ext->fetch(r, buffer.data());
        // Do something with each row.
    }
}, ptr->nrow(), num_threads);
```

`tatami_r::parallelize()` can also be used as a custom **tatami** parallelization scheme.
This requires a little bit of care as the `TATAMI_CUSTOM_PARALLEL` macro has to be set before any includes of **tatami**...
which itself is included by all **tatami_r** headers except for `parallelize.hpp`, so it's easy to get wrong!
We suggest using the following sequence of preprocessor statements before any **tatami**-containing source file.
(If using [**beachmat**](https://github.com/tatami-inc/beachmat)'s `Rtatami.h` header, all of these preprocessor statements have already been added for our convenience.)

```cpp
// Must be defined before including tatami_r/parallelize.hpp
#define TATAMI_R_PARALLELIZE_UNKNOWN 
#include "tatami_r/parallelize.hpp"

// Defining the tatami parallelization scheme, which must 
// occur before including tatami or tatami_r itself.
#define TATAMI_CUSTOM_PARALLEL tatami_r::parallelize

#include "tatami_r/tatami_r.hpp"
```

We can now use `tatami::parallelize()`, `TATAMI_CUSTOM_PARALLEL` and `tatami_r::parallelize()` interchangeably.

## Using the main thread executor

We can perform our own calls to the R API inside each worker by wrapping it in the **manticore** executor.
This will execute the user-provided function on the main thread before returning control to the worker.
Developers do not have to do the usual **manticore** dance of `initialize()`, `listen()`, `finish_thread()`, and so on; 
this is all handled by `tatami_r::parallelize()` itself, so only `run()` needs to be specified.

```cpp
auto& mexec = tatami_r::executor();

tatami_r::parallelize([&](int thread_id, int start, int len) -> void {
    mexec.run([&]() -> void {
        // Do something that touches the R API.
    });
}, ptr->nrow(), num_threads);
```

It is important to use the global executor provided by the `tatami_r::executor()` function, as this is the same as that used inside `tatami_r::parallelize()`.
Otherwise, if a different `manticore::Executor` instance is created, we will not be properly protected from simultaneous calls to the R API from different workers.

## Under the hood

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

Check out the implementation of `tatami_r::parallelize()` for more details.

## Dynamically loaded libraries

When working with multiple dynamically loaded libraries (e.g., in separate R packages), it is possible for `tatami_r::executor()` to return references to different instances.
(See [here](https://www.reddit.com/r/cpp_questions/comments/a5nhnm/why_is_a_static_function_variable_shared_between/) for a discussion on the relevant differences between Clang and GCC.)
This can be problematic if, e.g., a `tatami::Matrix` object is created in one library and then used in a parallel section in the other.

In such cases, we can force both libraries to use the same `manticore::Executor` instance.
We first obtain the address of the instance in one of the libraries, usually the one that is more upstream in the dependency chain:

```cpp
auto ptr = &(tatami::executor());
```

We then pass this address to the other library, using it to set that library's global instance.
This only needs to be done once, usually at library load time.

```cpp
tatami_r::set_executor(ptr);
```

After this is done, calls to `tatami_r::executor()` in the second library will use the instance of the first library.
This ensures that we maintain thread-safe calls to R across both libraries.

## Further comments

Construction of the `UnknownMatrix` must always be performed on the main thread.
This is because construction of the unknown fallback involves some calls into the R runtime; these are currently not protected from execution in worker contexts.
Similarly, any **Rcpp**-based allocations - even default construction of classes like `Rcpp::NumericVector` - should be done in the main thread, just in case.

`manticore::Executor::run()` is only required if the code can only be executed on the main thread.
For code that must be serial but does not need to be on the main thread, we can just use a standard synchronization primitives (e.g., `<mutex>`) inside `tatami_r::parallelize()` calls.
This works correctly as `tatami_r::parallelize()` uses `<thread>` under the hood.
