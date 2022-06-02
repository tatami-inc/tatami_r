# Convert R objects to `tatami` matrices

## Overview

**raticate** is an **Rcpp**-based header-only library for parsing R objects as `tatami::Matrix` instances.
This enables use of R matrix data in [**tatami**](https://github.com/LTLA/tatami)-compatible C++ code such as [**libscran**](https://github.com/LTLA/libscran).
Usage is as simple as:

```cpp
#include "raticate/raticate.hpp"

SEXP some_typical_rcpp_function(Rcpp::RObject x) {
    auto parsed = raticate::parse(x);

    if (parsed.matrix == std::nullptr) {
        // Do something if format of 'x' is not supported;
        // probably throw an error.
    } else {
        // Do stuff with the tatami::Matrix.
        parsed.matrix->nrow();
        auto first_row = parsed.matrix->row(0);
    }
}
```

And that's it, really.
If you want more details, you can check out the [reference documentation](https://ltla.github.io/raticate).

## Slightly more advanced usage

Currently `parse()` knows about the following matrix formats:

- ordinary logical, numeric or integer matrices.
- `dgCMatrix` or `lgCMatrix` objects from the **Matrix** package.
- `SparseArraySeed` objects from the **DelayedArray** package.
- `DelayedMatrix` objects wrapping any of the above, or containing the following delayed operations:
  - Subsetting
  - Modification of dimnames
  - Transposition

If `parse()` cannot interpret the format of `x`, the `.matrix` member will be set to a `nullptr`.
It is the caller's responsibility to handle this case.

The `parsed` output object will also contain an R list named `.contents`.
This holds references to the memory used by the `.matrix` member, preventing their premature garbage collection.
Callers should ensure that the `.matrix` member (or any copies) does not outlive the `.contents` R list.

## Deployment

**raticate** is intended to be compiled with other relevant C++ code inside an R package.
This is most easily done by modifying the package `DESCRIPTION` with:

```
LinkingTo: beachmat
```

which will automatically use the vendored copies of **raticate** (and **tatami**) inside the [**beachmat**](http://bioconductor.org/packages/beachmat) package.
Note that C++17 is required.

If **beachmat** cannot be used, then the R package developer will need to copy the **raticate** and **tatami** `include/` directories into the package's `inst/include`,
and then add a `Makevars` file like:

```
PKG_CPPFLAGS = -I../inst/include
```
