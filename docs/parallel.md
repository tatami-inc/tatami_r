# Handling parallelization

In general, **raticate** is thread-safe and can be used inside various parallelization constructs (e.g., OpenMP, `<thread>`).
This is true for all supported matrices listed in `parse()`, where the memory layout is well-defined.
At that point, **raticate** simply operates on raw pointers to the underlying data, without touching the R API at all.

The exception is when `parse()` is used with `allow_unknown = true`.
In this case, **raticate** needs to extract data from an unknown matrix fallback by calling `DelayedArray::extract_array()` via the R API.
The R API is strictly single-threaded, so we might naively think to lock all calls to R inside serial sections.
Unfortunately, this is not sufficient as we still get stack limit errors when R is called from a worker.
Some experimentation indicates that any call to R must only occur on the main thread.
I suspect that some other R-managed process (an event loop, perhaps?) is always running on the main thread;
calling R from the worker will fail to block this process and lead to parallel execution between the worker and the main threads.

The solution is to force all worker threads to request evaluation of R code from the main thread.
Specifically, when a worker is inside a parallel section, it will send a message to the main thread requesting evaluation of `extract_array()`.
The main thread will perform the evaluation, convert the results into something independent of R, and return those results to the worker.
This ensures that the R runtime is safely evaluated without requiring all **Rcpp** code to be lifted out of the `tatami::Matrix`'s getter methods
(which would otherwise require a complete redesign of all applications that consume `tatami::Matrix` objects),

While it is possible for users to implement their own parallelization scheme based on the logic above, it is probably easier to use the pre-packaged `raticate::parallelize()` method.
This uses C++11's parallelization facilities (e.g., `<thread>`, `<mutex>`) to implement the threading and message passing.
Note that the `RATICATE_PARALLELIZE_UNKNOWN` macro must be defined in order for the parallel-aware code to be defined;
otherwise **raticate** will just use the much simpler single-threaded code if parallelization is not required.

```cpp
// This macro must be defined before including raticate.
#define RATICATE_PARALLELIZE_UNKNOWN
#include "raticate/raticate.hpp"

tatami::Matrix<double, int>* mat; // created somewhere else.
size_t nthreads = 10;

raticate::parallelize<double, int>(
    mat->ncol(), 
    [&](size_t start, size_t end) -> void {
        for (size_t i = 0; i < end; ++i) {
            auto x = mat->column(i);
            // do something with the column.
        }
    }, 
    nthreads
);
```

The `raticate::parallelize()` function is most obviously applied in `tatami::apply()`.
This requires a bit of preprocessing trickery to ensure that the definitions are available at the right time.
It also, unfortunately, requires prior knowledge of the type of the `tatami::Matrix` - in this case, we are assuming that we are dealing with a `tatami::Matrix<double, int>`.
Also see `tests/src/bindings.cpp` for a working example.

```cpp
// First defining the macros before raticate or tatami are included.
template<class Function> 
void custom_run(size_t n, Function f);
#define TATAMI_CUSTOM_PARALLEL custom_run
#define RATICATE_PARALLELIZE_UNKNOWN

// Including raticate, which also includes tatami.
#include "raticate/raticate.hpp"

// Now actually defining the custom_run function.
template<class Function> 
void run(size_t n, Function f) {
    raticate::parallelize<double, int>(n, f, 3);
}
```

Note that `raticate::parse()` must always be called from the main thread.
This is because construction of the unknown fallback involves some calls into the R runtime;
these are currently not protected from execution in worker contexts.
