/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef PME_FFT_H__
#define PME_FFT_H__

#ifdef NAMD_FFTW
#include <rfftw.h>
#elif defined(NAMD_SGI_COMPLIB_FFT)
extern "C" {
#include <fft.h>
}
#endif

class PmeFFT {

public:
  PmeFFT(int K1, int K2, int K3);
  ~PmeFFT();
  
  // Actual physical dimensions of array
  void getdims(int *qdim2, int *qdim3) {*qdim2=dim2; *qdim3=dim3;}

  void forward(double *q_arr);
  void backward(double *q_arr);
  
private:
  int k1, k2, k3, dim2, dim3;
#if defined(NAMD_FFTW)
  rfftwnd_plan forward_plan;
  rfftwnd_plan backward_plan;
#elif defined(NAMD_SGI_COMPLIB_FFT)
  double *sgi_complib_fft_coeff;
#else
  int ntable;
  double *table;
  double *work;
#endif

};

#endif

