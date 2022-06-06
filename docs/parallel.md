# Handling parallelization

For the most part, **raticate** is thread-safe and can be used inside various parallelization constructs (e.g., OpenMP, `<thread>`).
The exception is when **raticate** needs to call the R runtime (via **Rcpp**), which is strictly single-threaded.
In such cases, **raticate** will attempt to enforce serial execution.

By default, it is assumed that parallelization is performed using OpenMP, so **raticate** wraps all **Rcpp** calls inside OpenMP `critical` sections.
These sections are named according to the `RATICATE_RCPP_CRITICAL_NAME` macro.
Advanced callers can ensure that **raticate** calls are executed in serial alongside other **Rcpp** calls by doing something like:

```cpp
// Define this before including raticate.
#define RATICATE_RCPP_CRITICAL_NAME rcpp

#include "raticate/raticate.hpp"

#pragma omp critical(rcpp)
{
    // Do something unrelated with Rcpp, which will be
    // thread-safe with any concurrent raticate calls.
}
```

For other parallelization schemes, users can define the `RATICATE_RCPP_PARALLEL_LOCK` macro.
This should be a template function that accepts a single lambda and executes it with no arguments.
For example, a caller parallelizing via `<thread>` might use the following to achieve  thread safety:

```cpp
#include <thread>

std::mutex global_rcpp_lock;

template<class Function>
void lock_for_rcpp(Function fun) {
    std::lock_guard<std::mutex> lock(global_rcpp_lock);
    fun();
}

// Define this before including raticate.
#define RATICATE_RCPP_PARALLEL_LOCK lock_for_rcpp
```
