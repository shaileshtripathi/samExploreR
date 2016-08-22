#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include <R.h>
#include "regrsub.h"

static const R_CMethodDef cMethods[]  = {
  {"R_readSummary_wrapper",(DL_FUNC)&"R_readSummary_wrapper",2},
  {NULL, NULL, 0}
};

void R_init_flagme(DllInfo *info){
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}
