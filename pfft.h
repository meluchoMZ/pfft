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
#include <thread>
#include <mutex>
#include <condition_variable>

#define ZERO_THRESHOLD 0.0000000001 // any number below 10^(-10) will be rounded to zero
#define FFT_TYPE_RECURSIVE 0x00
#define FFT_TYPE_ITERATIVE 0x00



// computes 1 dimension FFT using the classic recursive divide & conquer algorithm
// x is the input signal, y the output signal. Ideally, y will empty, cause its resized inside the function
// y length will be a power of two due to external padding (if needed)
void FFT_recursive(std::valarray<std::complex<double>> &x, std::valarray<std::complex<double>> &y);

// computes 1 dimension FFT in an iterative way
void FFT_iterative(std::valarray<std::complex<double>> &x, std::valarray<std::complex<double>> &y);

// computes 1 dimension IFFT using the 1D FFT classic recursive algorithm
// y is the input signal, x the output signal. Ideally, x will be empty, cause its resized inside the function
// y length must be a power of two 
// x external zero padding will be removed of the output sample set
void IFFT_recursive(std::valarray<std::complex<double>> &y, std::valarray<std::complex<double>> &x);

// computes 1 dimension IFFT in an iterative way
void IFFT_iterative(std::valarray<std::complex<double>> &y, std::valarray<std::complex<double>> &x);


// computes the 2D FFT. Only works for input lengths power of 2, cause is radix 2
void FFT2_recursive(std::valarray<std::complex<double>> &x, std::valarray<std::complex<double>> &y, unsigned long columns);
void FFT2_iterative(std::valarray<std::complex<double>> &x, std::valarray<std::complex<double>> &y, unsigned long columns);
// computes the 2D IFFT. Only works for input lengths power of 2, cause is radix 2
void IFFT2_recursive(std::valarray<std::complex<double>> &y, std::valarray<std::complex<double>> &x, unsigned long columns);
void IFFT2_iterative(std::valarray<std::complex<double>> &y, std::valarray<std::complex<double>> &x, unsigned long columns);


// parallel implementation of the FFT using shared memory
// thread_number must be a power of two, lower than min(nºcolumns, nºrows)
void PFFT2_recursive(std::valarray<std::complex<double>> &x, std::valarray<std::complex<double>> &y, unsigned long columns, int thread_number);
void PFFT2_iterative(std::valarray<std::complex<double>> &x, std::valarray<std::complex<double>> &y, unsigned long columns, int thread_number);

// parallel implementation of the IFFT using shared memory
// thread_number must be a power of two, lower than min(nºcolumns, nºrows)
void PIFFT2_recursive(std::valarray<std::complex<double>> &y, std::valarray<std::complex<double>> &x, unsigned long columns, int thread_number);
void PIFFT2_iterative(std::valarray<std::complex<double>> &y, std::valarray<std::complex<double>> &x, unsigned long columns, int thread_number);

#endif
