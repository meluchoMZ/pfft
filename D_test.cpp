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
	//__1D_FFT_recursive(data);
	FFT_recursive(data);
	std::cout << "FFT:" << std::endl;
	for (int p = 0; p < 17; p++) {
		std::cout << "FT[" << input_data[p] << "] = " << data[p] << std::endl;
	}
	std::valarray<std::complex<double>> fftd = data;
	__1D_IFFT_recursive(fftd);
	std::cout << "IFFT:" << std::endl;
	for (int p = 0; p < 17; p++) {
		std::cout << "IFFT[" << data[p] << "] = " << fftd[p] << std::endl;
	}
	return EXIT_SUCCESS;
}
