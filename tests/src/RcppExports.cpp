// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// parse
SEXP parse(Rcpp::RObject seed, double cache_size, bool require_min);
RcppExport SEXP _raticate_tests_parse(SEXP seedSEXP, SEXP cache_sizeSEXP, SEXP require_minSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< double >::type cache_size(cache_sizeSEXP);
    Rcpp::traits::input_parameter< bool >::type require_min(require_minSEXP);
    rcpp_result_gen = Rcpp::wrap(parse(seed, cache_size, require_min));
    return rcpp_result_gen;
END_RCPP
}
// num_rows
int num_rows(Rcpp::RObject parsed);
RcppExport SEXP _raticate_tests_num_rows(SEXP parsedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type parsed(parsedSEXP);
    rcpp_result_gen = Rcpp::wrap(num_rows(parsed));
    return rcpp_result_gen;
END_RCPP
}
// num_columns
int num_columns(Rcpp::RObject parsed);
RcppExport SEXP _raticate_tests_num_columns(SEXP parsedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type parsed(parsedSEXP);
    rcpp_result_gen = Rcpp::wrap(num_columns(parsed));
    return rcpp_result_gen;
END_RCPP
}
// prefer_rows
bool prefer_rows(Rcpp::RObject parsed);
RcppExport SEXP _raticate_tests_prefer_rows(SEXP parsedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type parsed(parsedSEXP);
    rcpp_result_gen = Rcpp::wrap(prefer_rows(parsed));
    return rcpp_result_gen;
END_RCPP
}
// is_sparse
bool is_sparse(Rcpp::RObject parsed);
RcppExport SEXP _raticate_tests_is_sparse(SEXP parsedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type parsed(parsedSEXP);
    rcpp_result_gen = Rcpp::wrap(is_sparse(parsed));
    return rcpp_result_gen;
END_RCPP
}
// myopic_dense_full
Rcpp::List myopic_dense_full(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx);
RcppExport SEXP _raticate_tests_myopic_dense_full(SEXP parsedSEXP, SEXP rowSEXP, SEXP idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type parsed(parsedSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type idx(idxSEXP);
    rcpp_result_gen = Rcpp::wrap(myopic_dense_full(parsed, row, idx));
    return rcpp_result_gen;
END_RCPP
}
// oracular_dense_full
Rcpp::List oracular_dense_full(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx);
RcppExport SEXP _raticate_tests_oracular_dense_full(SEXP parsedSEXP, SEXP rowSEXP, SEXP idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type parsed(parsedSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type idx(idxSEXP);
    rcpp_result_gen = Rcpp::wrap(oracular_dense_full(parsed, row, idx));
    return rcpp_result_gen;
END_RCPP
}
// myopic_dense_block
Rcpp::List myopic_dense_block(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, int first, int len);
RcppExport SEXP _raticate_tests_myopic_dense_block(SEXP parsedSEXP, SEXP rowSEXP, SEXP idxSEXP, SEXP firstSEXP, SEXP lenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type parsed(parsedSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< int >::type first(firstSEXP);
    Rcpp::traits::input_parameter< int >::type len(lenSEXP);
    rcpp_result_gen = Rcpp::wrap(myopic_dense_block(parsed, row, idx, first, len));
    return rcpp_result_gen;
END_RCPP
}
// oracular_dense_block
Rcpp::List oracular_dense_block(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, int first, int len);
RcppExport SEXP _raticate_tests_oracular_dense_block(SEXP parsedSEXP, SEXP rowSEXP, SEXP idxSEXP, SEXP firstSEXP, SEXP lenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type parsed(parsedSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< int >::type first(firstSEXP);
    Rcpp::traits::input_parameter< int >::type len(lenSEXP);
    rcpp_result_gen = Rcpp::wrap(oracular_dense_block(parsed, row, idx, first, len));
    return rcpp_result_gen;
END_RCPP
}
// myopic_dense_indexed
Rcpp::List myopic_dense_indexed(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, Rcpp::IntegerVector subset);
RcppExport SEXP _raticate_tests_myopic_dense_indexed(SEXP parsedSEXP, SEXP rowSEXP, SEXP idxSEXP, SEXP subsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type parsed(parsedSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type subset(subsetSEXP);
    rcpp_result_gen = Rcpp::wrap(myopic_dense_indexed(parsed, row, idx, subset));
    return rcpp_result_gen;
END_RCPP
}
// oracular_dense_indexed
Rcpp::List oracular_dense_indexed(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, Rcpp::IntegerVector subset);
RcppExport SEXP _raticate_tests_oracular_dense_indexed(SEXP parsedSEXP, SEXP rowSEXP, SEXP idxSEXP, SEXP subsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type parsed(parsedSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type subset(subsetSEXP);
    rcpp_result_gen = Rcpp::wrap(oracular_dense_indexed(parsed, row, idx, subset));
    return rcpp_result_gen;
END_RCPP
}
// myopic_sparse_full
Rcpp::List myopic_sparse_full(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, bool needs_value, bool needs_index);
RcppExport SEXP _raticate_tests_myopic_sparse_full(SEXP parsedSEXP, SEXP rowSEXP, SEXP idxSEXP, SEXP needs_valueSEXP, SEXP needs_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type parsed(parsedSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< bool >::type needs_value(needs_valueSEXP);
    Rcpp::traits::input_parameter< bool >::type needs_index(needs_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(myopic_sparse_full(parsed, row, idx, needs_value, needs_index));
    return rcpp_result_gen;
END_RCPP
}
// oracular_sparse_full
Rcpp::List oracular_sparse_full(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, bool needs_value, bool needs_index);
RcppExport SEXP _raticate_tests_oracular_sparse_full(SEXP parsedSEXP, SEXP rowSEXP, SEXP idxSEXP, SEXP needs_valueSEXP, SEXP needs_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type parsed(parsedSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< bool >::type needs_value(needs_valueSEXP);
    Rcpp::traits::input_parameter< bool >::type needs_index(needs_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(oracular_sparse_full(parsed, row, idx, needs_value, needs_index));
    return rcpp_result_gen;
END_RCPP
}
// myopic_sparse_block
Rcpp::List myopic_sparse_block(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, int first, int len, bool needs_value, bool needs_index);
RcppExport SEXP _raticate_tests_myopic_sparse_block(SEXP parsedSEXP, SEXP rowSEXP, SEXP idxSEXP, SEXP firstSEXP, SEXP lenSEXP, SEXP needs_valueSEXP, SEXP needs_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type parsed(parsedSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< int >::type first(firstSEXP);
    Rcpp::traits::input_parameter< int >::type len(lenSEXP);
    Rcpp::traits::input_parameter< bool >::type needs_value(needs_valueSEXP);
    Rcpp::traits::input_parameter< bool >::type needs_index(needs_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(myopic_sparse_block(parsed, row, idx, first, len, needs_value, needs_index));
    return rcpp_result_gen;
END_RCPP
}
// oracular_sparse_block
Rcpp::List oracular_sparse_block(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, int first, int len, bool needs_value, bool needs_index);
RcppExport SEXP _raticate_tests_oracular_sparse_block(SEXP parsedSEXP, SEXP rowSEXP, SEXP idxSEXP, SEXP firstSEXP, SEXP lenSEXP, SEXP needs_valueSEXP, SEXP needs_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type parsed(parsedSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< int >::type first(firstSEXP);
    Rcpp::traits::input_parameter< int >::type len(lenSEXP);
    Rcpp::traits::input_parameter< bool >::type needs_value(needs_valueSEXP);
    Rcpp::traits::input_parameter< bool >::type needs_index(needs_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(oracular_sparse_block(parsed, row, idx, first, len, needs_value, needs_index));
    return rcpp_result_gen;
END_RCPP
}
// myopic_sparse_indexed
Rcpp::List myopic_sparse_indexed(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, Rcpp::IntegerVector subset, bool needs_value, bool needs_index);
RcppExport SEXP _raticate_tests_myopic_sparse_indexed(SEXP parsedSEXP, SEXP rowSEXP, SEXP idxSEXP, SEXP subsetSEXP, SEXP needs_valueSEXP, SEXP needs_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type parsed(parsedSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type subset(subsetSEXP);
    Rcpp::traits::input_parameter< bool >::type needs_value(needs_valueSEXP);
    Rcpp::traits::input_parameter< bool >::type needs_index(needs_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(myopic_sparse_indexed(parsed, row, idx, subset, needs_value, needs_index));
    return rcpp_result_gen;
END_RCPP
}
// oracular_sparse_indexed
Rcpp::List oracular_sparse_indexed(Rcpp::RObject parsed, bool row, Rcpp::IntegerVector idx, Rcpp::IntegerVector subset, bool needs_value, bool needs_index);
RcppExport SEXP _raticate_tests_oracular_sparse_indexed(SEXP parsedSEXP, SEXP rowSEXP, SEXP idxSEXP, SEXP subsetSEXP, SEXP needs_valueSEXP, SEXP needs_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type parsed(parsedSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type subset(subsetSEXP);
    Rcpp::traits::input_parameter< bool >::type needs_value(needs_valueSEXP);
    Rcpp::traits::input_parameter< bool >::type needs_index(needs_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(oracular_sparse_indexed(parsed, row, idx, subset, needs_value, needs_index));
    return rcpp_result_gen;
END_RCPP
}
// myopic_dense_sums
Rcpp::NumericVector myopic_dense_sums(Rcpp::RObject parsed, bool row, int num_threads);
RcppExport SEXP _raticate_tests_myopic_dense_sums(SEXP parsedSEXP, SEXP rowSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type parsed(parsedSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(myopic_dense_sums(parsed, row, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// oracular_dense_sums
Rcpp::NumericVector oracular_dense_sums(Rcpp::RObject parsed, bool row, int num_threads);
RcppExport SEXP _raticate_tests_oracular_dense_sums(SEXP parsedSEXP, SEXP rowSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type parsed(parsedSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(oracular_dense_sums(parsed, row, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// myopic_sparse_sums
Rcpp::NumericVector myopic_sparse_sums(Rcpp::RObject parsed, bool row, int num_threads);
RcppExport SEXP _raticate_tests_myopic_sparse_sums(SEXP parsedSEXP, SEXP rowSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type parsed(parsedSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(myopic_sparse_sums(parsed, row, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// oracular_sparse_sums
Rcpp::NumericVector oracular_sparse_sums(Rcpp::RObject parsed, bool row, int num_threads);
RcppExport SEXP _raticate_tests_oracular_sparse_sums(SEXP parsedSEXP, SEXP rowSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type parsed(parsedSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(oracular_sparse_sums(parsed, row, num_threads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_raticate_tests_parse", (DL_FUNC) &_raticate_tests_parse, 3},
    {"_raticate_tests_num_rows", (DL_FUNC) &_raticate_tests_num_rows, 1},
    {"_raticate_tests_num_columns", (DL_FUNC) &_raticate_tests_num_columns, 1},
    {"_raticate_tests_prefer_rows", (DL_FUNC) &_raticate_tests_prefer_rows, 1},
    {"_raticate_tests_is_sparse", (DL_FUNC) &_raticate_tests_is_sparse, 1},
    {"_raticate_tests_myopic_dense_full", (DL_FUNC) &_raticate_tests_myopic_dense_full, 3},
    {"_raticate_tests_oracular_dense_full", (DL_FUNC) &_raticate_tests_oracular_dense_full, 3},
    {"_raticate_tests_myopic_dense_block", (DL_FUNC) &_raticate_tests_myopic_dense_block, 5},
    {"_raticate_tests_oracular_dense_block", (DL_FUNC) &_raticate_tests_oracular_dense_block, 5},
    {"_raticate_tests_myopic_dense_indexed", (DL_FUNC) &_raticate_tests_myopic_dense_indexed, 4},
    {"_raticate_tests_oracular_dense_indexed", (DL_FUNC) &_raticate_tests_oracular_dense_indexed, 4},
    {"_raticate_tests_myopic_sparse_full", (DL_FUNC) &_raticate_tests_myopic_sparse_full, 5},
    {"_raticate_tests_oracular_sparse_full", (DL_FUNC) &_raticate_tests_oracular_sparse_full, 5},
    {"_raticate_tests_myopic_sparse_block", (DL_FUNC) &_raticate_tests_myopic_sparse_block, 7},
    {"_raticate_tests_oracular_sparse_block", (DL_FUNC) &_raticate_tests_oracular_sparse_block, 7},
    {"_raticate_tests_myopic_sparse_indexed", (DL_FUNC) &_raticate_tests_myopic_sparse_indexed, 6},
    {"_raticate_tests_oracular_sparse_indexed", (DL_FUNC) &_raticate_tests_oracular_sparse_indexed, 6},
    {"_raticate_tests_myopic_dense_sums", (DL_FUNC) &_raticate_tests_myopic_dense_sums, 3},
    {"_raticate_tests_oracular_dense_sums", (DL_FUNC) &_raticate_tests_oracular_dense_sums, 3},
    {"_raticate_tests_myopic_sparse_sums", (DL_FUNC) &_raticate_tests_myopic_sparse_sums, 3},
    {"_raticate_tests_oracular_sparse_sums", (DL_FUNC) &_raticate_tests_oracular_sparse_sums, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_raticate_tests(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
