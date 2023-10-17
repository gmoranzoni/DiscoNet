#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

void mcode_complex(int*, int*, float*, float*, int*, int*, int*);

static const R_CMethodDef CEntries[] = {
  {"mcode_complex", (DL_FUNC) &mcode_complex, 7},
  {NULL, NULL, 0}
};

void R_init_DiscoNet(DllInfo *dll) {
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
