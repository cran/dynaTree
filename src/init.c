#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void alc_R(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void alcX_R(void *, void *, void *, void *, void *, void *);
extern void classprobs_R(void *, void *, void *, void *, void *, void *, void *);
extern void coef_R(void *, void *, void *, void *, void *);
extern void copy_cloud_R(void *, void *);
extern void delete_cloud_R(void *);
extern void delete_clouds_R();
extern void dynaTree_R(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void entropyX_R(void *, void *, void *);
extern void ieci_R(void *, void *, void *, void *, void *, void *, void *, void *);
extern void intervals_R(void *, void *, void *, void *, void *);
extern void predclass_R(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void predict_R(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void qEI_R(void *, void *, void *, void *, void *, void *, void *);
extern void qEntropy_R(void *, void *, void *, void *, void *, void *);
extern void rejuvenate_R(void *, void *, void *, void *, void *);
extern void relevance_R(void *, void *, void *, void *, void *, void *);
extern void retire_R(void *, void *, void *, void *, void *, void *, void *);
extern void sameleaf_R(void *, void *, void *, void *);
extern void sens_R(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void treestats_R(void *, void *, void *, void *, void *);
extern void update_R(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void varproptotal_R(void *, void *);
extern void varpropuse_R(void *, void *);

static const R_CMethodDef CEntries[] = {
    {"alc_R",           (DL_FUNC) &alc_R,           10},
    {"alcX_R",          (DL_FUNC) &alcX_R,           6},
    {"classprobs_R",    (DL_FUNC) &classprobs_R,     7},
    {"coef_R",          (DL_FUNC) &coef_R,           5},
    {"copy_cloud_R",    (DL_FUNC) &copy_cloud_R,     2},
    {"delete_cloud_R",  (DL_FUNC) &delete_cloud_R,   1},
    {"delete_clouds_R", (DL_FUNC) &delete_clouds_R,  0},
    {"dynaTree_R",      (DL_FUNC) &dynaTree_R,      14},
    {"entropyX_R",      (DL_FUNC) &entropyX_R,       3},
    {"ieci_R",          (DL_FUNC) &ieci_R,           8},
    {"intervals_R",     (DL_FUNC) &intervals_R,      5},
    {"predclass_R",     (DL_FUNC) &predclass_R,      9},
    {"predict_R",       (DL_FUNC) &predict_R,       16},
    {"qEI_R",           (DL_FUNC) &qEI_R,            7},
    {"qEntropy_R",      (DL_FUNC) &qEntropy_R,       6},
    {"rejuvenate_R",    (DL_FUNC) &rejuvenate_R,     5},
    {"relevance_R",     (DL_FUNC) &relevance_R,      6},
    {"retire_R",        (DL_FUNC) &retire_R,         7},
    {"sameleaf_R",      (DL_FUNC) &sameleaf_R,       4},
    {"sens_R",          (DL_FUNC) &sens_R,          21},
    {"treestats_R",     (DL_FUNC) &treestats_R,      5},
    {"update_R",        (DL_FUNC) &update_R,        10},
    {"varproptotal_R",  (DL_FUNC) &varproptotal_R,   2},
    {"varpropuse_R",    (DL_FUNC) &varpropuse_R,     2},
    {NULL, NULL, 0}
};

void R_init_dynaTree(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
