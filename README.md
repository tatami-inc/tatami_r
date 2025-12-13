# Read R objects via tatami 

![Unit tests](https://github.com/tatami-inc/tatami_r/actions/workflows/run-tests.yaml/badge.svg)
![Documentation](https://github.com/tatami-inc/tatami_r/actions/workflows/doxygenate.yaml/badge.svg)

## Overview

**tatami_r** is an header-only library for reading abstract R matrices in [**tatami**](https://github.com/tatami-inc/tatami).
This allows **tatami**-based C++ functions to accept and operate on any matrix-like R object containing numeric data.
Usage is as simple as:

```cpp
#include "tatami_r/tatami_r.hpp"

SEXP some_typical_rcpp_function(Rcpp::RObject x) {
    auto ptr = std::make_shared<tatami_r::UnknownMatrix<double, int> >(x);

    // Do stuff with the tatami::Matrix.
    ptr->nrow();
    auto row_extractor = ptr->dense_row();
    auto first_row = row_extractor->fetch(0);

    // Return something.
    return R_NilValue;
}
```

For more details, check out the [reference documentation](https://tatami-inc.github.io/tatami_r).

## Implementation

**tatami_r** assumes that the hosting R instance has loaded the [**DelayedArray**](https://bioconductor.org/packages/DelayedArray) package.
The `UnknownMatrix` will then use the `extract_array()` and `extract_sparse_array()` R functions to retrieve data from the abstract R matrix.
Note that this involves calling into R from C++, so high performance should not be expected here.
Rather, the purpose of **tatami_r** is to ensure that **tatami**-based functions keep working when a native representation cannot be found for a particular matrix-like object.

Most R package developers will not need to use **tatami_r** directly.
Rather, they should use the `initializeCpp()` function from the [**beachmat**](https://bioconductor.org/packages/beachmat) package to map an arbitrary matrix to its appropriate representation.
When such mappings exist, this allows the C++ code to operate without calling back into R for maximum efficiency.
If no mapping is known, **beachmat** will gracefully fall back to an `UnknownMatrix` to keep things running.

# Enabling parallelization

Given a `tatami_r::UnknownMatrix` or a `tatami::Matrix*` that might refer to one, we can easily parallelize operations with the `tatami_r::parallelize()` function.
This accepts a lambda/functor with the thread ID and the range of jobs (in the example below, rows) to be processed.

```cpp
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

Any calls to the `extract_*_array()` R functions are made thread-safe by the [**manticore**](https://github.com/tatami-inc/manticore) library.
Developers can also access the **manticore** executor to safely perform their own R API calls from each thread.

```cpp
auto& mexec = tatami_r::executor();

tatami_r::parallelize([&](size_t thread_id, int start, int len) -> void {
    mexec.run([&]() -> void {
        // Do something that touches the R API.
    });
}, ptr->nrow(), num_threads);
```

Check out the [comments about safe parallelization](docs/parallel.md) for more gory details.

## Deployment

**tatami_r** is intended to be compiled with other relevant C++ code inside an R package using [**Rcpp**](https://www.rcpp.org/).
This is most easily done by modifying the package `DESCRIPTION` with:

```
LinkingTo: beachmat, assorthead, Rcpp
```

which will automatically use the vendored copies of **tatami_r** (and **tatami**) inside the [**assorthead**](http://bioconductor.org/packages/assorthead) package,
along with some of pre-configured macro definitions for safe parallelization in [**beachmat**](https://bioconductor.org/packages/beachmat)'s `Rtatami.h` header.
Note that C++17 is required.

If **assorthead** or **beachmat** cannot be used, developers should ensure that the contents of the `include/` directories
(as well as all dependencies listed in [`extern/CMakeLists.txt`](extern/CMakeLists.txt))
are available during package build, and then add a `Makevars` file like:

```
PKG_CPPFLAGS = -Isome/path/to/headers
```
