// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppGSL.h>
#include <Rcpp.h>

using namespace Rcpp;

// simulateFixedTreeCodeNew
int simulateFixedTreeCodeNew(int nObs, int traverseFrequency, std::string indir, std::string outdir, SEXP seed);
RcppExport SEXP _diffpop_simulateFixedTreeCodeNew(SEXP nObsSEXP, SEXP traverseFrequencySEXP, SEXP indirSEXP, SEXP outdirSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nObs(nObsSEXP);
    Rcpp::traits::input_parameter< int >::type traverseFrequency(traverseFrequencySEXP);
    Rcpp::traits::input_parameter< std::string >::type indir(indirSEXP);
    Rcpp::traits::input_parameter< std::string >::type outdir(outdirSEXP);
    Rcpp::traits::input_parameter< SEXP >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(simulateFixedTreeCodeNew(nObs, traverseFrequency, indir, outdir, seed));
    return rcpp_result_gen;
END_RCPP
}
// simulateTreeCodeNew
int simulateTreeCodeNew(int nObs, int traverseFrequency, std::string indir, std::string outdir, SEXP seed);
RcppExport SEXP _diffpop_simulateTreeCodeNew(SEXP nObsSEXP, SEXP traverseFrequencySEXP, SEXP indirSEXP, SEXP outdirSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nObs(nObsSEXP);
    Rcpp::traits::input_parameter< int >::type traverseFrequency(traverseFrequencySEXP);
    Rcpp::traits::input_parameter< std::string >::type indir(indirSEXP);
    Rcpp::traits::input_parameter< std::string >::type outdir(outdirSEXP);
    Rcpp::traits::input_parameter< SEXP >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(simulateTreeCodeNew(nObs, traverseFrequency, indir, outdir, seed));
    return rcpp_result_gen;
END_RCPP
}
// testcode
int testcode();
RcppExport SEXP _diffpop_testcode() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(testcode());
    return rcpp_result_gen;
END_RCPP
}
// testcode2
int testcode2(int l);
RcppExport SEXP _diffpop_testcode2(SEXP lSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    rcpp_result_gen = Rcpp::wrap(testcode2(l));
    return rcpp_result_gen;
END_RCPP
}
// test3
Rcpp::NumericVector test3(int size, int numReps);
RcppExport SEXP _diffpop_test3(SEXP sizeSEXP, SEXP numRepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< int >::type numReps(numRepsSEXP);
    rcpp_result_gen = Rcpp::wrap(test3(size, numReps));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_diffpop_simulateFixedTreeCodeNew", (DL_FUNC) &_diffpop_simulateFixedTreeCodeNew, 5},
    {"_diffpop_simulateTreeCodeNew", (DL_FUNC) &_diffpop_simulateTreeCodeNew, 5},
    {"_diffpop_testcode", (DL_FUNC) &_diffpop_testcode, 0},
    {"_diffpop_testcode2", (DL_FUNC) &_diffpop_testcode2, 1},
    {"_diffpop_test3", (DL_FUNC) &_diffpop_test3, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_diffpop(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
