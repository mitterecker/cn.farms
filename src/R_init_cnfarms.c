#include "cn.farms.h"
#include <R_ext/Rdynload.h>

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

/* laplace.cpp */
	CALLMETHOD_DEF(momentsGauss, 10),

/* sparse_farms.c */
	CALLMETHOD_DEF(sparseFarmsC, 4),
	CALLMETHOD_DEF(getL, 1),
	CALLMETHOD_DEF(getEss, 1),
	CALLMETHOD_DEF(getLap, 1),
	CALLMETHOD_DEF(deinit, 1),

	{NULL, NULL, 0}
};


void R_init_cnfarms(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	return;
}
