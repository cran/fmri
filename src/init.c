#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
#include <Rinternals.h> // for SEXP
#include <R_ext/RS.h>

void F77_NAME(ccluster)(int* x, int* n1, int* n2, int* n3, int* z);
void F77_NAME(extrpatt)(double* beta, int* voxel, int* n1, int* n2, int* n3,
                      int* nb, int* sl, int* nsl, double* pattern, int* nvox);
void F77_NAME(getslpv)(double* stat, int* n, double* p, double* kv, int* nsim,
                       double* pval);
void F77_NAME(lconnect)(int* segm, int* n1, int* n2, int* n3, int* i1, int* i2, int* i3,
                        int* ind1, int* ind2, int* ind3, int* checked, int* mask);
void F77_NAME(segm3d)(double* y, double* res, double* si2, int* pos,
                      int* wlse, int* n1, int* n2, int* n3, int* nt, double* df,
                      double* hakt, double* lambda, double* theta,
                      double* bi, double* thn, double* lwght, double* wght,
                      double* swres, double* pval, int* segm, double* delta,
                      double* thresh, double* fov, double* vq,
                      double* vest0i, double* varest, int* restrict);
void F77_NAME(sincfilter)(double* t, int* nt, double* x, int* nx, double* ft, int* wr);
void F77_NAME(slicetim)(double* x, int* nt, int* n1, int* n2, int* n3, double* y,
                        double* t, int* sliceord);
void F77_NAME(slight)(double* stat, int* mask, int* n1, int* n2, int* n3,
                      int* slght, int* nsl, double* slstat);
void F77_NAME(smooth3d)(double* y, double* si2, int* mask, int* wlse, int* n1,
                        int* n2, int* n3, int* dv, double* hakt, double* thn,
                        int* kern, double* lwght, double* wght, double* swjy);
void F77_NAME(thcorr)(double* w, int* n1, int* n2, int* n3, double* corr,
                      int* l1, int* l2, int* l3);

static R_NativePrimitiveArgType ccluster_t[]={INTSXP, INTSXP, INTSXP,
    INTSXP, INTSXP};
static R_NativePrimitiveArgType extrpatt_t[]={REALSXP, INTSXP, INTSXP,
    INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType getslpv_t[]={REALSXP, INTSXP, REALSXP,
    REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType lconnect_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
    INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType segm3d_t[]={REALSXP, REALSXP, REALSXP, INTSXP,
    INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType sincfilter_t[]={REALSXP, INTSXP, REALSXP,
    INTSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType slicetim_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
    INTSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType slight_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
    INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType thcorr_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
    REALSXP, INTSXP, INTSXP, INTSXP};


static const R_FortranMethodDef fmethods[] = {
            {"ccluster", (DL_FUNC) &ccluster_ ,5 , ccluster_t},
            {"extrpatt", (DL_FUNC) &extrpatt_ ,10, extrpatt_t},
            {"getslpv", (DL_FUNC) &getslpv_ ,6, getslpv_t},
            {"lconnect", (DL_FUNC) &lconnect_ ,12, lconnect_t},
            {"segm3d", (DL_FUNC) &segm3d_ ,27, segm3d_t},
            {"sincfilter", (DL_FUNC) &sincfilter_ ,6, sincfilter_t},
            {"slicetim", (DL_FUNC) &slicetim_ ,8, slicetim_t},
            {"slight", (DL_FUNC) &slight_ ,8, slight_t},
            {"thcorr", (DL_FUNC) &thcorr_ ,8, thcorr_t},
           {NULL, NULL, 0,NULL}
};

void R_init_fmri(DllInfo *dll)
         {
             R_registerRoutines(dll, NULL, NULL, fmethods , NULL);
             R_useDynamicSymbols(dll,FALSE);
             R_forceSymbols(dll,TRUE);
         }
