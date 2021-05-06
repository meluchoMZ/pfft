/**
  Computes FFT and 2D FFT through serveral methods
  Digital information processing, Computer engineering
Author: Miguel Blanco Godón
*/

#ifndef __PFFT_MELUCHO_MZ
#define __PFFT_MELUCHO_MZ

#include <valarray>
#include <complex>
#include <cmath>
#include <iostream>

// computes 1 dimension FFT using the classic recursive divide & conquer algorithm
void FFT_recursive(std::valarray<std::complex<double>> &x);
void __1D_FFT_recursive(std::valarray<std::complex<double>> &x);

// computes 1 dimension IFFT using the 1D FFT classic recursive algorithm
void __1D_IFFT_recursive(std::valarray<std::complex<double>> &X);

#endif
