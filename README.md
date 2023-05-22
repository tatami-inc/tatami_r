# Read R objects via tatami 

## Overview

**tatami_r** is an header-only library for reading matrix-like R objects in [**tatami**](https://github.com/tatami-inc/tatami).
Usage is as simple as:

```cpp
#include "tatami_r/tatami_r.hpp"

SEXP some_typical_rcpp_function(Rcpp::RObject x) {
    auto ptr = tatami_r::UnknownMatrix(x);

    // Do stuff with the tatami::Matrix.
    ptr->nrow();
    auto row_extractor = ptr->dense_row();
    auto first_row = row_extractor->fetch(0);
}
```

And that's it, really.
If you want more details, you can check out the [reference documentation](https://tatami-inc.github.io/tatami_r).
Also check out the [comments about safe parallelization](docs/parallel.md) when dealing with `tatami::Matrix` pointers that might contain `tatami_r::UnknownMatrix` objects.

## Deployment

**tatami_r** is intended to be compiled with other relevant C++ code inside an R package.
This is most easily done by modifying the package `DESCRIPTION` with:

```
LinkingTo: beachmat
```

which will automatically use the vendored copies of **tatami_r** (and **tatami**) inside the [**beachmat**](http://bioconductor.org/packages/beachmat) package.
Note that C++17 is required.

If **beachmat** cannot be used, then the R package developer will need to copy the **tatami_r** and **tatami** `include/` directories into the package's `inst/include`,
and then add a `Makevars` file like:

```
PKG_CPPFLAGS = -I../inst/include
```
