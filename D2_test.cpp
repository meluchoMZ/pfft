/**
  Test for the 2D FFT and IFFT
  Digital information processing, Computer engineering
Author: Miguel Blanco Godón
*/

#include <iostream>
#include <complex>
#include "pfft.h"


int main(void)
{
	//std::complex<double> input_data[] = {1, 5, 10, 15, 2, 4, 6, 8};
	std::complex<double> input_data[] = {2, -4, 5, 6, 12, -5, -3, 2, 5, 7, -1, 3, -8, -9, 6, 4};
	std::valarray<std::complex<double>> data(input_data, 16);
	std::valarray<std::complex<double>> fftd, recv;

	std::cout << "Input matrix:" << std::endl;
	for (int k = 0; k < 8; k++) {
		for (int kk = 0; kk < 2; kk++) {
			std::cout << data[k*2 + kk] << " ";
		}
		std::cout << std::endl;
	}

	FFT2(data, fftd, 2, FFT_TYPE_RECURSIVE);

	std::cout << std::endl << "FFT matrix:" << std::endl;
	for (int k = 0; k < 8; k++) {
		for (int kk = 0; kk < 2; kk++) {
			std::cout << fftd[k*2+kk] << " ";
		}
		std::cout << std::endl;
	}

	IFFT2(fftd, recv, 2, FFT_TYPE_ITERATIVE);

	std::cout << std::endl << "IFFT matrix:" << std::endl;
	for (int k = 0; k < 8; k++) {
		for (int kk = 0; kk < 2; kk++) {
			std::cout << recv[k*2+kk] << " ";
		}
		std::cout << std::endl;
	}

	std::valarray<std::complex<double>> copy;
	copy = recv;
	copy.resize(0);
	
	std::cout << "copy size: " << copy.size() << std::endl;
	std::cout << "recv (orig) size: " << recv.size() << std::endl;
	//recv[1:10];
	return EXIT_SUCCESS;
}
