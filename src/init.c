#include <R.h>
#include <Rinternals.h>
#include "Rdefines.h"
#include "R_ext/Rdynload.h"

void nodeExisted(double *nodeXlist,double *nodeYlist,int *n,double *X,double *Y, int *tag);
void edgelength(double *nodeXlist,double *nodeYlist,int *n, double *edgelength, int *longlat);
void gc_el(double *lon1, double *lon2, double *lat1, double *lat2, double *gel);
void dist_p2e(double *nodeXlist,double *nodeYlist,int *n,int *idx,double *X,double *Y, double *dist, double *fx, double *fy);
void footxy(double *nodeXlist,double *nodeYlist,double *X,double *Y, double *dist, double *fx, double *fy);
int maxindex(double *list, int n);
int minindex(double *list, int n);
void addDegree(int *NIDs, int *n,int *nid, int *DL);
void minusDegree(int *NIDs, int *n,int *nid, int *DL);

static const R_CMethodDef CEntries[] = {
    {"nodeExisted",              (DL_FUNC) &nodeExisted,              6},
    {"edgelength",          (DL_FUNC) &edgelength,          5},
    {"gc_el",       (DL_FUNC) &gc_el,       5},
    {"dist_p2e",       (DL_FUNC) &dist_p2e,       9},
    {"footxy",       (DL_FUNC) &footxy,       7},
    {"maxindex",      (DL_FUNC) &maxindex,      2},
    {"minindex",       (DL_FUNC) &minindex,       2},
    {"addDegree",            (DL_FUNC) &addDegree,            4},
    {"minusDegree", (DL_FUNC) &minusDegree, 4},
    {NULL, NULL, 0}
};

void R_init_shp2graph(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
