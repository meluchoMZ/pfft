/**
  Test FFT algorithms
  Digital information processing, Computer engineering
Author: Miguel Blanco Godón
*/

// C++ headers
#include <iostream>
#include <cstring>
#include <vector>
#include <complex>
// C headers
#include <dirent.h>
#include <errno.h>
#include <unistd.h>


int main(void)
{
	DIR *d; 
	struct dirent *dir;
	FILE *f;
	char c = 0x0;
	int pixel;
	std::complex<double> complex_pixel (0, 0);
	char test_dir[] = "./testbench";
	std::vector<unsigned char> pixels;
	std::vector<std::complex<double>> fft_vector;

	if (!(d = opendir(test_dir))) {
		std::cout << "Can't open '" << test_dir  << "' directory: " << strerror(errno) << std::endl;
		return EXIT_FAILURE;
	}
	while ((dir = readdir(d)) != NULL)
	{
		chdir(test_dir);
		if (strcmp(dir->d_name, ".") != 0 && strcmp(dir->d_name, "..") != 0) {
			std::cout << "Executing testbench " << c << ": " << dir->d_name << std::endl;
			if (!(f = fopen(dir->d_name, "r"))) {
				std::cout << "Can't open file '" << dir->d_name << "': " << strerror(errno) << std::endl;
				break;
			}
			while (fscanf(f, "%d,", &pixel) != EOF)
			{
				pixels.push_back((unsigned char) pixel);
				complex_pixel = std::complex<double> (pixel, 0);
				fft_vector.push_back(complex_pixel);
			}
			// FFT and IFFT application to the vector
			// value checking
			for (unsigned int i = 0; i < pixels.size(); i++) {
				if (pixels[i] != fft_vector[i].real()) {
					std::cout << "Computation differs in position " << i << ": " << fft_vector[i].real() << " != " << pixels[i] << std::endl;
				}
			}
			c++;
			fclose(f);
		}
		chdir("..");
	}
	closedir(d);
	std::cout << "Vector size: " << pixels.size() << std::endl;
	for (c = 0x0; c < 10; c++) {
		std::cout << (unsigned int) pixels[c] << " " <<  fft_vector[c] << std::endl;
	}
	pixels.clear();
	fft_vector.clear();
	std::cout << "Vector size after clearing: " << pixels.size() << std::endl;
	return EXIT_SUCCESS;
}
