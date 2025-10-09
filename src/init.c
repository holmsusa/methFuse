#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Declare your C functions exposed to R
extern SEXP sort_tree_R(SEXP);
extern SEXP fuse_cluster_R(SEXP, SEXP, SEXP, SEXP);
extern SEXP cuttree_R(SEXP, SEXP);
extern SEXP bino_div_R(SEXP, SEXP);
extern SEXP bino_xlogy_R(SEXP, SEXP);

// Register the functions for R's use
static const R_CallMethodDef CallEntries[] = {
  {"sort_tree_R", (DL_FUNC) &sort_tree_R, 1},
  {"fuse_cluster_R", (DL_FUNC) &fuse_cluster_R, 4},
  {"cuttree_R", (DL_FUNC) &cuttree_R, 2},
  {"bino_div_R", (DL_FUNC) &bino_div_R, 2},
  {"bino_xlogy_R", (DL_FUNC) &bino_xlogy_R, 2},
  {NULL, NULL, 0}
};

void R_init_fuseR(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE); // disables dynamic symbol lookup
}
