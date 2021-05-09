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

#define ZERO_THRESHOLD 0.0000000001 // any number below 10^(-10) will be rounded to zero

// computes 1 dimension FFT using the classic recursive divide & conquer algorithm
// x is the input signal, y the output signal. Ideally, y will empty, cause its resized inside the function
// y length will be a power of two due to external padding (if needed)
void FFT_recursive(std::valarray<std::complex<double>> &x, std::valarray<std::complex<double>> &y);

void FFT_iterative(std::valarray<std::complex<double>> &x, std::valarray<std::complex<double>> &y);

// computes 1 dimension IFFT using the 1D FFT classic recursive algorithm
// y is the input signal, x the output signal. Ideally, x will be empty, cause its resized inside the function
// y length must be a power of two 
// x external zero padding will be removed of the output sample set
void bit_reverse_copy(std::valarray<std::complex<double>> &original, std::valarray<std::complex<double>> &copy);

void IFFT_recursive(std::valarray<std::complex<double>> &y, std::valarray<std::complex<double>> &x);

void IFFT_iterative(std::valarray<std::complex<double>> &y, std::valarray<std::complex<double>> &x);

#endif
