/**
  Computes the FFT and 2D FFT 
  Digital information processing, Computer engineering
Author: Miguel Blanco Godón
*/

#include "pfft.h"

void get_new_size(unsigned long prev, unsigned long *next)
{
	*next = 2;
	while (*next < prev)
	{
		*next *= 2;
	}
}

void FFT_recursive(std::valarray<std::complex<double>> &x)
{
	unsigned long samples = x.size();
	unsigned long new_size;
	get_new_size(samples, &new_size);
	std::valarray<std::complex<double>> resized(0,0);
	resized.resize(new_size);
	resized[std::slice(0, samples, 1)] = x;
	std::cout << "RESIZED: " << std::endl;
	for (unsigned long k = 0; k < new_size; k++) {
		std::cout << k << ": " << resized[k] << std::endl;
	}
	__1D_FFT_recursive(resized);
	std::cout << "RESIZED FFT: " << std::endl;
	for (unsigned long k = 0; k < new_size; k++) {
		std::cout << k << ": " << resized[k] << std::endl;
	}
	x = resized[std::slice(0, samples, 1)];
}

void __1D_FFT_recursive(std::valarray<std::complex<double>> &x)
{
	unsigned long samples = x.size();
	if (samples <= 1) {
		return;
	} else {
		std::valarray<std::complex<double>> x_k_even = x[std::slice(0, samples/2, 2)];
		std::valarray<std::complex<double>> x_k_odd = x[std::slice(1, samples/2, 2)];
		__1D_FFT_recursive(x_k_even);
		__1D_FFT_recursive(x_k_odd);

		for (unsigned long k = 0; k < samples/2; k++) {
			std::complex<double> q = std::complex<double>(cos(-2*M_PI*k/samples), sin(-2*M_PI*k/samples)) * x_k_odd[k];
			x[k] = x_k_even[k]+q;
			x[k+samples/2] = x_k_even[k]-q;
		}
	}
}

void IFFT_recursive(std::valarray<std::complex<double>> &x)
{
}


void __1D_IFFT_recursive(std::valarray<std::complex<double>> &X)
{
	unsigned long samples = X.size();
	for (unsigned long k = 0; k < samples; k++) {
		X[k] = std::complex<double>(std::real(X[k]), -1*std::imag(X[k]));
	}
	__1D_FFT_recursive(X);
	for (unsigned long k = 0; k < samples; k++) {
		X[k] = std::complex<double>(std::real(X[k]), -1*std::imag(X[k]));
	}
	X /= samples;
}
