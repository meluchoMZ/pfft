/**
  Computes the FFT and 2D FFT 
  Digital information processing, Computer engineering
Author: Miguel Blanco Godón
*/

#include "pfft.h"


void get_new_size(unsigned long prev, unsigned long *next);
void __1D_FFT_recursive(std::valarray<std::complex<double>> &x);
void FFT_recursive(std::valarray<std::complex<double>> &x, std::valarray<std::complex<double>> &y);
void __1D_IFFT_recursive(std::valarray<std::complex<double>> &X);
void IFFT_recursive(std::valarray<std::complex<double>> &y, std::valarray<std::complex<double>> &x);
void bit_reverse_copy(std::valarray<std::complex<double>> &original, std::valarray<std::complex<double>> &copy);
void FFT_iterative(std::valarray<std::complex<double>> &x, std::valarray<std::complex<double>> &y);
void IFFT_iterative(std::valarray<std::complex<double>> &y, std::valarray<std::complex<double>> &x);
void __1D_FFT_iterative(std::valarray<std::complex<double>> &x);
void __1D_IFFT_iterative(std::valarray<std::complex<double>> &X);

// calculates the nearest power of two of prev and steores it in next
void get_new_size(unsigned long prev, unsigned long *next)
{
	*next = 2;
	while (*next < prev)
	{
		*next *= 2;
	}
}

/* Recursive versions */
/* 1D */ 
void FFT_recursive(std::valarray<std::complex<double>> &x, std::valarray<std::complex<double>> &y)
{
	unsigned long samples = x.size();
	unsigned long new_size;
	if (samples == 0) {
		return;
	}
	// replicates x on y, were length(y) is a power of two
	// applies zero padding if necessary
	get_new_size(samples, &new_size);
	y.resize(new_size);
	y[std::slice(0, samples, 1)] = x;
	// computes the fft
	__1D_FFT_recursive(y);
}

void __1D_FFT_recursive(std::valarray<std::complex<double>> &x)
{
	unsigned long samples = x.size();
	std::complex<double> q;

	if (samples <= 1) {
		return;
	} else {
		// divide the array in two parts:
		// 1. containts the even index elements of x
		// 2. containts the odd index elements of x
		std::valarray<std::complex<double>> x_k_even = x[std::slice(0, samples/2, 2)];
		std::valarray<std::complex<double>> x_k_odd = x[std::slice(1, samples/2, 2)];
		// computes de fft on each subvector
		__1D_FFT_recursive(x_k_even);
		__1D_FFT_recursive(x_k_odd);
		// gathers the computations and collects them back on the array
		for (unsigned long k = 0; k < samples/2; k++) {
			q = std::complex<double>(cos(-2*M_PI*k/samples), sin(-2*M_PI*k/samples)) * x_k_odd[k];
			x[k] = x_k_even[k]+q;
			x[k+samples/2] = x_k_even[k]-q;
		}
	}
}

void IFFT_recursive(std::valarray<std::complex<double>> &y, std::valarray<std::complex<double>> &x)
{
	unsigned long samples = y.size(), k = samples-1;
	std::valarray<std::complex<double>> aux;
	if (samples == 0) {
		return;
	}
	aux.resize(samples);
	// replicates y in an auxiliar vector
	aux[std::slice(0, samples, 1)] = y;
	// computes the ifft
	__1D_IFFT_recursive(aux);

	// removes zero padding and copies data to x
	while (abs(std::real(aux[k])) < ZERO_THRESHOLD && abs(std::imag(aux[k])) < ZERO_THRESHOLD)
	{
		k--;
	}
	x.resize(k+1);
	x[std::slice(0, k+1, 1)] = aux[std::slice(0, k+1, 1)];
}


void __1D_IFFT_recursive(std::valarray<std::complex<double>> &X)
{
	unsigned long samples = X.size();
	// computes complex conjugate of the input
	for (unsigned long k = 0; k < samples; k++) {
		X[k] = std::complex<double>(std::real(X[k]), -1*std::imag(X[k]));
	}
	// computes the fft
	__1D_FFT_recursive(X);
	// computes again the complex conjugate of the output of the fft
	for (unsigned long k = 0; k < samples; k++) {
		X[k] = std::complex<double>(std::real(X[k]), -1*std::imag(X[k]));
	}
	// normalizes the output
	X /= samples;
}


/* Iterative versions */

// copies original onto copy by bitwise reversing the index values
void bit_reverse_copy(std::valarray<std::complex<double>> &original, std::valarray<std::complex<double>> &copy)
{
	unsigned long s = original.size(), rev_k, rev_k_size, rev_num;
	copy.resize(s);
	// computes the length of the bit words
	if (s%2 == 0) {
		rev_k_size = log2(s);
	} else {
		rev_k_size = log2(s)+1;
	}
	for (unsigned long k = 0; k < s; k++) {
		rev_k = 0x00;
		rev_num = k;
		// this works like a circular shift register. v.gr
		// original: 0111
		// step 1: 1011
		// step 2: 1101
		// step 3: 1110
		for (unsigned long kk = 0; kk < rev_k_size; kk++) {
			rev_k = rev_k << 0x01;
			if ((rev_num & 0x01) == 0x01) {
				rev_k = rev_k ^ 0x01;
			}
			rev_num = rev_num >> 0x01;
		}
		// copies the values in the permutaded order
		copy[rev_k] = original[k];
	}
}


void FFT_iterative(std::valarray<std::complex<double>> &x, std::valarray<std::complex<double>> &y)
{
	unsigned long samples, new_size;
	std::valarray<std::complex<double>> aux;
	// replicates x and applies zero padding
	samples = x.size();
	if (samples == 0) {
		return;
	}
	get_new_size(samples, &new_size);
	aux.resize(new_size);
	aux[std::slice(0, samples, 1)] = x;
	bit_reverse_copy(aux, y);
	// computes the fft
	__1D_FFT_iterative(y);
}


void IFFT_iterative(std::valarray<std::complex<double>> &y, std::valarray<std::complex<double>> &x)
{
	unsigned long samples, k;
	std::valarray<std::complex<double>> aux, aux2;
	// replicates x and applies zero padding
	samples = y.size();
	if (samples == 0) {
		return;
	}
	k = samples -1;
	aux.resize(samples);
	aux[std::slice(0, samples, 1)] = y;
	bit_reverse_copy(aux, aux2);

	// computes the fft
	__1D_IFFT_iterative(aux2);

	// removes zero padding and copies data to x
	while (abs(std::real(aux2[k])) < ZERO_THRESHOLD && abs(std::imag(aux2[k])) < ZERO_THRESHOLD)
	{
		k--;
	}
	x.resize(k+1);
	x[std::slice(0, k+1, 1)] = aux2[std::slice(0, k+1, 1)];
}


void __1D_FFT_iterative(std::valarray<std::complex<double>> &y)
{
	unsigned long samples, lsize, m;
	std::complex<double> w, w_m, t, u;
	samples = y.size();
	lsize = log2(samples);
	for (unsigned long s = 1; s <= lsize; s++) {
		m = 1 << s;
		w_m = std::complex<double>(cos(-2*M_PI/m), sin(-2*M_PI/m));
		for (unsigned long k = 0; k < samples; k += m) {
			w = 1;
			// in situ computation of the butterflies
			for (unsigned long j = 0; j < (m >> 1); j++) {
				t = w * y[k + j + (m >> 1)];
				u = y[k + j];
				y[k + j] = u + t;
				y[k + j + (m >> 1)] = u - t;
				w *= w_m;
			}
		}
	}
}

void __1D_IFFT_iterative(std::valarray<std::complex<double>> &X)
{
	unsigned long samples = X.size();
	// computes complex conjugate of the input
	for (unsigned long k = 0; k < samples; k++) {
		X[k] = std::complex<double>(std::real(X[k]), -1*std::imag(X[k]));
	}
	// computes the fft
	__1D_FFT_iterative(X);
	// computes again the complex conjugate of the output of the fft
	for (unsigned long k = 0; k < samples; k++) {
		X[k] = std::complex<double>(std::real(X[k]), -1*std::imag(X[k]));
	}
	// normalizes the output
	X /= samples;
}
