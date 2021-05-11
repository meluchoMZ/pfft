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

	//FFT2_recursive(data, fftd, 2);
	//FFT2_iterative(data, fftd, 2);
	//PFFT2_recursive(data, fftd, 2, 2);
	PFFT2_iterative(data, fftd, 2, 2);

	std::cout << std::endl << "FFT matrix:" << std::endl;
	for (int k = 0; k < 8; k++) {
		for (int kk = 0; kk < 2; kk++) {
			std::cout << fftd[k*2+kk] << " ";
		}
		std::cout << std::endl;
	}

	//IFFT2_recursive(fftd, recv, 2);
	IFFT2_iterative(fftd, recv, 2);

	std::cout << std::endl << "IFFT matrix:" << std::endl;
	for (int k = 0; k < 8; k++) {
		for (int kk = 0; kk < 2; kk++) {
			std::cout << recv[k*2+kk] << " ";
		}
		std::cout << std::endl;
	}
	
	return EXIT_SUCCESS;
}
