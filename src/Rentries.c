#include <R.h>
#include <Rmath.h>

double F77_SUB(fpchisq)(double *x, double *df, int *low_tail, int *give_log) {
   return pchisq(*x,*df,*low_tail,*give_log);
}

