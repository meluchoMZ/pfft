/**
  Test FFT algorithms
  Digital information processing, Computer engineering
Author: Miguel Blanco Godón
*/

#include <iostream>
#include <cstring>
#include <complex>
#include <vector>
#include <valarray>

#include <dirent.h>
#include <errno.h>
#include <unistd.h>
#include <sys/time.h>

#include "pfft.h"

#define TEST0 "./testbench/test_0.dat"
#define TEST1 "./testbench/test_1.dat"
#define TEST2 "./testbench/test_2.dat"
#define TEST3 "./testbench/test_3.dat"
#define TEST4 "./testbench/test_4.dat"

int main(void)
{
	FILE *f;
	int d;
	std::vector<std::complex<double>> input_vector;
	std::valarray<std::complex<double>> *data;
	std::valarray<std::complex<double>> fftd;
	struct timeval start, finish;
	std::string test_files[] = {TEST0, TEST1, TEST2, TEST3, TEST4}; 
	int columns[] = {256, 512, 1024, 4096, 8192};


	for (int k = 0; k < 5; k++) {
		if ((f = fopen(test_files[k].c_str(), "r")) == NULL) {
			std::cout << "ERROR: cannot open test file '" << test_files[k].c_str() << "': " << strerror(errno) << "." << std::endl;
		}
	
		while (fscanf(f, "%d\n", &d) != EOF)
		{
			input_vector.push_back(std::complex<double> (d, 0));
		}
		
		data = new std::valarray<std::complex<double>> (input_vector.data(), input_vector.size());
		
		gettimeofday(&start, NULL);
		PFFT2_recursive(*data, fftd, columns[k], 4);
		gettimeofday(&finish, NULL);
		std::cout << "Execution time of 2D PFFT (recursive) on test " << k << "(" << columns[k] << "x" << columns[k] << "): " << (double) ((finish.tv_sec-start.tv_sec)*1000+(finish.tv_usec-start.tv_usec)/1000)/1000 << " seconds." << std::endl;

		gettimeofday(&start, NULL);
		PFFT2_iterative(*data, fftd, columns[k], 4);
		gettimeofday(&finish, NULL);
		std::cout << "Execution time of 2D PFFT (iterative) on test " << k << "(" << columns[k] << "x" << columns[k] << "): " << (double) ((finish.tv_sec-start.tv_sec)*1000+(finish.tv_usec-start.tv_usec)/1000)/1000 << " seconds." << std::endl;
		delete data;

		fclose(f); f = NULL;
		input_vector.clear();
	}

	return EXIT_SUCCESS;
}
