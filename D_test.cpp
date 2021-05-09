/**
  1D FFT and IFFT test
  Digital information processing, computer engineering
Author: Miguel Blanco Godón
*/


#include <iostream>
#include <complex>
#include "pfft.h"

int main(void)
{
	std::complex<double> input_data[] = {1, 2, 3, 4, 5, 6, 7, 8, -1, -2, -3, -4, -5, -6, -7, -8, -9};
	std::valarray<std::complex<double>> data(input_data, 17);
	std::valarray<std::complex<double>> fftd, recovered, fftdi, recoveredi;
	std::cout << "Input data:" << std::endl;
	for (int p = 0; p < (int)data.size(); p++) {
		std::cout << data[p] << std::endl;
	}
	FFT_recursive(data, fftd);
	std::cout << std::endl << "FFT recursive:" << std::endl;
	for (int p = 0; p < (int)fftd.size(); p++) {
		std::cout << fftd[p] << std::endl;
	}
	IFFT_recursive(fftd, recovered);
	std::cout << std::endl << "IFFT:" << std::endl;
	for (int p = 0; p < (int)recovered.size(); p++) {
		std::cout << recovered[p] << std::endl;
	}

	FFT_iterative(data, fftdi);
	std::cout << std::endl << "FFT iterative:" << std::endl;
	for (int p = 0; p < (int)fftdi.size(); p++) {
		std::cout << fftdi[p] << std::endl;
	}
	IFFT_iterative(fftdi, recoveredi);
	std::cout << std::endl << "IFFT iterative:" << std::endl;
	for (int p = 0; p < (int)recoveredi.size(); p++) {
		std::cout << recoveredi[p] << std::endl;
	}
	return EXIT_SUCCESS;
}
