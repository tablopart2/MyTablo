#ifndef QRNG_SOBOL_H
#define QRNG_SOBOL_H

void qrng_sobol(int dim,int length, double* seq); 
//dim: dimension
//length: length of each dimension
//seq: resulting sequence, size of seq should be dim*lenght*sizeof(double)
#endif