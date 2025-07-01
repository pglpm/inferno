#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "initseq.h"

static R_CMethodDef cMethods[] = {
    {NULL, NULL, 0, NULL}
};

static R_CallMethodDef callMethods[]  = {
    {"initseq", (DL_FUNC) &initseq, 1},
    {NULL, NULL, 0}
};

void attribute_visible R_init_mcmc(DllInfo *info)
{
    R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}
