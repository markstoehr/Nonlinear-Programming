#include <stdio.h>

// multiply a sparse symmetric Hessian
void symm_matvec_mult(int *Hridx, int *Hcidx, double *Hval, double *inVec,double *outVec,int nIdx)
{
  
  int ridx,cidx;
  for (int i=0; i<nIdx; i++) {
    ridx = Hridx[i];
    cidx = Hcidx[i];
    if (ridx != cidx) {
      outVec[ridx]+=Hval[i]*inVec[cidx];
      outVec[cidx]+=Hval[i]*inVec[ridx];
    } else 
      outVec[ridx]+=Hval[i]*inVec[cidx];
  }
  
}
