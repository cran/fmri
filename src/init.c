#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
#include <Rinternals.h> // for SEXP
#include <R_ext/RS.h>

void F77_NAME(ccluster)(int* x, int* n1, int* n2, int* n3, int* z);
void F77_NAME(chaws2)(double* y, double* si2, int* mask, int* wlse, int* n1,
                      int* n2, int* n3, double* hakt, double* lambda,
                      double* theta, double* bi, double* thn, int* kern, int* skern,
                      double* spmin, double* spmax, double* lwght, double* wght);
void F77_NAME(chawsv)(double* y, double* res, double* si2, int* mask, int* wlse,
                      int* n1, int* n2, int* n3, int* n4, double* hakt,
                      double* lambda, double* theta, double* bi, double* resnew,
                      double* thn, int* kern, int* skern, double* spmin,
                      double* spmax, double* lwght, double* wght, double* resi);
void F77_NAME(extrpatt)(double* beta, int* voxel, int* n1, int* n2, int* n3,
                      int* nb, int* sl, int* nsl, double* pattern, int* nvox);
void F77_NAME(gethani)(double* x, double* y, int* kern, double* value,
                       double* wght, double* eps, double* bw);
void F77_NAME(getslpv)(double* stat, int* n, double* p, double* kv, int* nsim,
                       double* pval);
void F77_NAME(getvofh)(double* bw, int* kern, double* wght, double* vol);
void F77_NAME(ihaws2)(double* y, double* si2, int* mask, int* wlse, int* n1,
                      int* n2, int* n3, int* dv, double* hakt, double* lambda,
                      double* theta, int* ncores, double* bi, double* thn,
                      int* kern, int* skern, double* spmin, double* spmax,
                      double* lwght, double* wght, double* swjy);
void F77_NAME(imcorr)(double* res, int* mask, int* n1, int* n2, int* n3, int* nv,
                      double* scorr, int* l1, int* l2, int* l3);
void F77_NAME(ivar)(double* res, double* resscale, int* mask, int* n1, int* n2,
                    int* n3, int* nv, double* var);
void F77_NAME(lconnect)(int* segm, int* n1, int* n2, int* n3, int* i1, int* i2, int* i3,
                        int* ind1, int* ind2, int* ind3, int* checked, int* mask);
void F77_NAME(mcorr)(double* res, int* mask, int* n1, int* n2, int* n3, int* nv,
                     double* scorr, int* l1, int* l2, int* l3);
void F77_NAME(mean3d)(double* res, int* n1, int* n2, int* n3, int* nv, double* mres);
void F77_NAME(segm3d)(double* y, double* res, double* si2, int* mask,
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
void F77_NAME(sofw3df)(double* bw, int* kern, double* wght, double* fw);
void F77_NAME(thcorr)(double* w, int* n1, int* n2, int* n3, double* corr,
                      int* l1, int* l2, int* l3);

static R_NativePrimitiveArgType ccluster_t[]={INTSXP, INTSXP, INTSXP,
    INTSXP, INTSXP};
static R_NativePrimitiveArgType chaws2_t[]={REALSXP, REALSXP, LGLSXP,
    LGLSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP,
    REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType chawsv_t[]={REALSXP, REALSXP, REALSXP,
    LGLSXP, LGLSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType extrpatt_t[]={REALSXP, LGLSXP, INTSXP,
    INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType gethani_t[]={REALSXP, REALSXP, INTSXP,
    REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType getslpv_t[]={REALSXP, INTSXP, REALSXP,
    REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType getvofh_t[]={REALSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType ihaws2_t[]={REALSXP, REALSXP, LGLSXP, LGLSXP,
    INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP,
    REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType imcorr_t[]={REALSXP, LGLSXP, INTSXP, INTSXP,
    INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType ivar_t[]={REALSXP, REALSXP, LGLSXP, INTSXP,
    INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType lconnect_t[]={LGLSXP, INTSXP, INTSXP, INTSXP,
    INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, LGLSXP, LGLSXP};
static R_NativePrimitiveArgType mcorr_t[]={REALSXP, LGLSXP, INTSXP, INTSXP,
    INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType mean3d_t[]={REALSXP, INTSXP, INTSXP,
    INTSXP, INTSXP, REALSXP };
static R_NativePrimitiveArgType segm3d_t[]={REALSXP, REALSXP, REALSXP, LGLSXP,
    LGLSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, LGLSXP};
static R_NativePrimitiveArgType sincfilter_t[]={REALSXP, INTSXP, REALSXP,
    INTSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType slicetim_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
    INTSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType slight_t[]={REALSXP, LGLSXP, INTSXP, INTSXP,
    INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType smooth3d_t[]={REALSXP, REALSXP, LGLSXP, LGLSXP,
    INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP,
    REALSXP};
static R_NativePrimitiveArgType sofw3df_t[]={REALSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType thcorr_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
    REALSXP, INTSXP, INTSXP, INTSXP};


static const R_FortranMethodDef fmethods[] = {
            {"ccluster", (DL_FUNC) &ccluster_ ,5 , ccluster_t},
            {"chaws2", (DL_FUNC) &chaws2_ ,18, chaws2_t},
            {"chawsv", (DL_FUNC) &chawsv_ ,22, chawsv_t},
            {"extrpatt", (DL_FUNC) &extrpatt_ ,10, extrpatt_t},
            {"gethani", (DL_FUNC) &gethani_ ,7, gethani_t},
            {"getslpv", (DL_FUNC) &getslpv_ ,6, getslpv_t},
            {"getvofh", (DL_FUNC) &getvofh_ ,4, getvofh_t},
            {"ihaws2", (DL_FUNC) &ihaws2_ ,21, ihaws2_t},
            {"imcorr", (DL_FUNC) &imcorr_ ,10, imcorr_t},
            {"ivar", (DL_FUNC) &ivar_ ,8, ivar_t},
            {"lconnect", (DL_FUNC) &lconnect_ ,12, lconnect_t},
            {"mcorr", (DL_FUNC) &mcorr_ ,10, mcorr_t},
            {"mean3d", (DL_FUNC) &mean3d_ ,6, mean3d_t},
            {"segm3d", (DL_FUNC) &segm3d_ ,27, segm3d_t},
            {"sincfilter", (DL_FUNC) &sincfilter_ ,6, sincfilter_t},
            {"slicetim", (DL_FUNC) &slicetim_ ,8, slicetim_t},
            {"slight", (DL_FUNC) &slight_ ,8, slight_t},
            {"smooth3d", (DL_FUNC) &smooth3d_ ,14, smooth3d_t},
            {"sofw3df", (DL_FUNC) &sofw3df_ ,4, sofw3df_t},
            {"thcorr", (DL_FUNC) &thcorr_ ,8, thcorr_t},
           {NULL, NULL, 0,NULL}
};

void R_init_fmri(DllInfo *dll)
         {
             R_registerRoutines(dll, NULL, NULL, fmethods , NULL);
             R_useDynamicSymbols(dll,FALSE);
             R_forceSymbols(dll,TRUE);
         }
