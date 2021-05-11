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
	*next = 0x02;
	while (*next < prev)
	{
		*next <<= 1;
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
	if (samples % 2 != 0) {
		get_new_size(samples, &new_size);
		y.resize(new_size);
		y[std::slice(0, samples, 1)] = x;
	} else {
		y = x;
	}
	// computes the fft
	__1D_FFT_recursive(y);
}

void __1D_FFT_recursive(std::valarray<std::complex<double>> &x)
{
	unsigned long samples = x.size();
	std::complex<double> q;
	register unsigned long msamples = samples >> 1;

	if (samples <= 1) {
		return;
	} else {
		// divide the array in two parts:
		// 1. containts the even index elements of x
		// 2. containts the odd index elements of x
		std::valarray<std::complex<double>> x_k_even = x[std::slice(0, msamples, 2)];
		std::valarray<std::complex<double>> x_k_odd = x[std::slice(1, msamples, 2)];
		// computes de fft on each subvector
		__1D_FFT_recursive(x_k_even);
		__1D_FFT_recursive(x_k_odd);
		// gathers the computations and collects them back on the array
		for (register unsigned long k = 0; k < msamples; k++) {
			//q = std::complex<double>(cos(-2*M_PI*k/samples), sin(-2*M_PI*k/samples)) * x_k_odd[k];
			q = std::polar(1.0, -M_PI*(k<<1)/samples) * x_k_odd[k];
			x[k] = x_k_even[k]+q;
			x[k+msamples] = x_k_even[k]-q;
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
	double l2 = log2(s);
	//copy.resize(s);
	copy = original;
	//std::cout << "vector size: " << s << std::endl;
	// computes the length of the bit words
	if ((int) l2 == (double) l2) {
		rev_k_size = l2;
	} else {
		rev_k_size = (int) l2+1;
	}
	for (register unsigned long k = 0; k < s; k++) {
		rev_k = 0x00;
		rev_num = k;
		// this works like a circular shift register. v.gr
		// original: 0111
		// step 1: 1011
		// step 2: 1101
		// step 3: 1110
		for (register unsigned long kk = 0; kk < rev_k_size; kk++) {
			rev_k = rev_k << 0x01;
			if ((rev_num & 0x01) == 0x01) {
				rev_k = rev_k ^ 0x01;
			}
			rev_num = rev_num >> 0x01;
		}
		// copies the values in the permutaded order
		//std::cout << "original k: " << k << ";   permutaded k: " << rev_k << std::endl;
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
	std::valarray<std::complex<double>> aux;
	// replicates x and applies zero padding
	samples = y.size();
	if (samples == 0) {
		return;
	}
	k = samples -1;
	aux = y;

	bit_reverse_copy(y, aux);

	// computes the fft
	__1D_IFFT_iterative(aux);

	// removes zero padding and copies data to x
	while (abs(std::real(aux[k])) < ZERO_THRESHOLD && abs(std::imag(aux[k])) < ZERO_THRESHOLD)
	{
		k--;
	}
	x.resize(k+1);
	x[std::slice(0, k+1, 1)] = aux[std::slice(0, k+1, 1)];
}


void __1D_FFT_iterative(std::valarray<std::complex<double>> &y)
{
	unsigned long samples, lsize, m;
	std::complex<double> w, w_m, t, u;
	samples = y.size();
	lsize = log2(samples);
	for (unsigned long s = 1; s <= lsize; s++) {
		m = 1 << s;
		//w_m = std::complex<double>(cos(-2*M_PI/m), sin(-2*M_PI/m));
		w_m = std::polar(1.0, -M_PI/(m >> 1));
		for (register unsigned long k = 0; k < samples; k += m) {
			w = 1;
			// in situ computation of the butterflies
			for (register unsigned long j = 0; j < (m >> 1); j++) {
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
	register unsigned long samples = X.size();
	// computes complex conjugate of the input
	for (register unsigned long k = 0; k < samples; k++) {
		X[k] = std::complex<double>(std::real(X[k]), -1*std::imag(X[k]));
	}
	// computes the fft
	__1D_FFT_iterative(X);
	// computes again the complex conjugate of the output of the fft
	for (register unsigned long k = 0; k < samples; k++) {
		X[k] = std::complex<double>(std::real(X[k]), -1*std::imag(X[k]));
	}
	// normalizes the output
	X /= samples;
}



// 2D FFT's

void FFT2_recursive(std::valarray<std::complex<double>> &x, std::valarray<std::complex<double>> &y, unsigned long columns)
{
	unsigned long samples;
	register unsigned long rows;
	std::valarray<std::complex<double>> aux;
	samples = x.size();
	if (samples == 0) {
		return;
	}
	y = x;
	rows = samples/columns;
	aux.resize(rows);
	// computes fft on columns
	for (register unsigned long k = 0; k < columns; k++) {
		aux[std::slice(0, rows, 1)] = x[std::slice(k, rows, columns)];
		__1D_FFT_recursive(aux);
		y[std::slice(k, rows, columns)] = aux[std::slice(0, rows, 1)];
	}
	aux.resize(columns);
	// computes fft on rows
	for (register unsigned long k = 0; k < rows; k++) {
		aux[std::slice(0, columns, 1)] = y[std::slice(k*columns, columns, 1)];
		__1D_FFT_recursive(aux);
		y[std::slice(k*columns, columns, 1)] = aux[std::slice(0, columns, 1)];
	}
}


void IFFT2_recursive(std::valarray<std::complex<double>> &y, std::valarray<std::complex<double>> &x, unsigned long columns)
{
	unsigned long samples, rows;
	std::valarray<std::complex<double>> aux;
	samples = y.size();
	if (samples == 0) {
		return;
	}
	x.resize(samples);
	rows = samples/columns;
	aux.resize(rows);
	// computes ifft on columns
	for (unsigned long k = 0; k < columns; k++) {
		aux[std::slice(0, rows, 1)] = y[std::slice(k, rows, columns)];
		__1D_IFFT_recursive(aux);
		x[std::slice(k, rows, columns)] = aux[std::slice(0, rows, 1)];
	}
	aux.resize(columns);
	// computes ifft on rows
	for (unsigned long k = 0; k < rows; k++) {
		aux[std::slice(0, columns, 1)] = x[std::slice(k*columns, columns, 1)];
		__1D_IFFT_recursive(aux);
		x[std::slice(k*columns, columns, 1)] = aux[std::slice(0, columns, 1)];
	}
}



void FFT2_iterative(std::valarray<std::complex<double>> &x, std::valarray<std::complex<double>> &y, unsigned long columns)
{
	unsigned long samples;
	register unsigned long rows;
	std::valarray<std::complex<double>> aux, aux2;
	samples = x.size();
	if (samples == 0) {
		return;
	}
	y = x;
	rows = samples/columns;
	aux.resize(rows);
	// computes fft on columns
	for (register unsigned long k = 0; k < columns; k++) {
		aux[std::slice(0, rows, 1)] = x[std::slice(k, rows, columns)];
		bit_reverse_copy(aux, aux2);
		__1D_FFT_iterative(aux2);
		y[std::slice(k, rows, columns)] = aux2[std::slice(0, rows, 1)];
	}
	aux.resize(columns);
	// computes fft on rows
	for (register unsigned long k = 0; k < rows; k++) {
		aux[std::slice(0, columns, 1)] = y[std::slice(k*columns, columns, 1)];
		bit_reverse_copy(aux,aux2);
		__1D_FFT_iterative(aux2);
		y[std::slice(k*columns, columns, 1)] = aux2[std::slice(0, columns, 1)];
	}
}


void IFFT2_iterative(std::valarray<std::complex<double>> &x, std::valarray<std::complex<double>> &y, unsigned long columns)
{
	unsigned long samples;
	register unsigned long rows;
	std::valarray<std::complex<double>> aux, aux2;
	samples = x.size();
	if (samples == 0) {
		return;
	}
	y = x;
	rows = samples/columns;
	aux.resize(rows);
	// computes fft on columns
	for (register unsigned long k = 0; k < columns; k++) {
		aux[std::slice(0, rows, 1)] = x[std::slice(k, rows, columns)];
		bit_reverse_copy(aux, aux2);
		__1D_IFFT_iterative(aux2);
		y[std::slice(k, rows, columns)] = aux2[std::slice(0, rows, 1)];
	}
	aux.resize(columns);
	// computes fft on rows
	for (register unsigned long k = 0; k < rows; k++) {
		aux[std::slice(0, columns, 1)] = y[std::slice(k*columns, columns, 1)];
		bit_reverse_copy(aux,aux2);
		__1D_IFFT_iterative(aux2);
		y[std::slice(k*columns, columns, 1)] = aux2[std::slice(0, columns, 1)];
	}
}



// Parallel 2D FFT's

void recursive_column_wrapper(std::valarray<std::complex<double>> &x, std::valarray<std::complex<double>> &y, unsigned long alpha, unsigned long beta, unsigned long gamma, int who, std::mutex &m)
{
	std::valarray<std::complex<double>> aux;
	//m.lock();
	//m.unlock();
	//m.try_lock();
	//m.unlock();
	aux.resize(alpha);
	for (register unsigned long k = 0; k < gamma; k++) {
		aux[std::slice(0, alpha, 1)] = x[std::slice(k+who*gamma, alpha, beta)];
		__1D_FFT_recursive(aux);
		
		while(1)
		{
			if (m.try_lock()) {
				y[std::slice(k+who*gamma, alpha, beta)] = aux[std::slice(0, alpha, 1)];
				m.unlock();
				//cv.notify_all();
				break;
			} 
			//cv.wait(m);
		}
		//y[std::slice(k+who*gamma, alpha, beta)] = aux[std::slice(0, alpha, 1)];
		//*/
	}
	return;
}

void recursive_row_wrapper(std::valarray<std::complex<double>> &x, std::valarray<std::complex<double>> &y, unsigned long alpha, unsigned long beta, unsigned long gamma, int who, std::mutex &m)
{
	std::valarray<std::complex<double>> aux;
	aux.resize(beta);
	alpha += 0;
	for (unsigned long k = 0; k < gamma; k++) {
		//m.lock();
		//std::cout << "I am " << who << ". std::slice parameters: " << beta*(k+who*gamma) << " " << beta << " " << 1 << std::endl;
		//m.unlock();
		aux[std::slice(0, beta, 1)] = x[std::slice(beta*(k+who*gamma), beta, 1)];
		__1D_FFT_recursive(aux);
		while (1)
		{
			if (m.try_lock()) {
				y[std::slice(beta*(k+who*gamma), beta, 1)] = aux[std::slice(0, beta, 1)];
				m.unlock();
				break;
			}
		}
	}
}

void PFFT2_recursive(std::valarray<std::complex<double>> &x, std::valarray<std::complex<double>> &y, unsigned long columns, int thread_number)
{
	unsigned long samples, rows;
	std::valarray<std::complex<double>> aux;
	std::thread threads[thread_number];
	std::mutex m;

	samples = x.size();
	if (samples == 0) {
		return;
	}
	aux.resize(samples);
	y.resize(samples);
	rows = samples/columns;
	
	for (int k = 0; k < thread_number; k++) {
		threads[k] = std::thread(recursive_column_wrapper, std::ref(x), std::ref(aux), rows, columns, columns/thread_number, k, std::ref(m));
	}	
	for (int k = 0; k < thread_number; k++) {
		threads[k].join();
	}

	for (int k = 0; k < thread_number; k++) {
		threads[k] = std::thread(recursive_row_wrapper, std::ref(aux), std::ref(y), rows, columns, rows/thread_number, k, std::ref(m));
	}	
	for (int k = 0; k < thread_number; k++) {
		threads[k].join();
	}
}





void iterative_column_wrapper(std::valarray<std::complex<double>> &x, std::valarray<std::complex<double>> &y, unsigned long alpha, unsigned long beta, unsigned long gamma, int who, std::mutex &m)
{
	std::valarray<std::complex<double>> aux, aux2;
	aux.resize(alpha);
	for (register unsigned long k = 0; k < gamma; k++) {
		aux[std::slice(0, alpha, 1)] = x[std::slice(k+who*gamma, alpha, beta)];
		bit_reverse_copy(aux, aux2);
		__1D_FFT_recursive(aux2);
		
		while(1)
		{
			if (m.try_lock()) {
				y[std::slice(k+who*gamma, alpha, beta)] = aux2[std::slice(0, alpha, 1)];
				m.unlock();
				//cv.notify_all();
				break;
			} 
			//cv.wait(m);
		}
		//y[std::slice(k+who*gamma, alpha, beta)] = aux[std::slice(0, alpha, 1)];
		//*/
	}
	return;
}

void iterative_row_wrapper(std::valarray<std::complex<double>> &x, std::valarray<std::complex<double>> &y, unsigned long alpha, unsigned long beta, unsigned long gamma, int who, std::mutex &m)
{
	std::valarray<std::complex<double>> aux, aux2;
	aux.resize(beta);
	alpha += 0;
	for (unsigned long k = 0; k < gamma; k++) {
		//m.lock();
		//std::cout << "I am " << who << ". std::slice parameters: " << beta*(k+who*gamma) << " " << beta << " " << 1 << std::endl;
		//m.unlock();
		aux[std::slice(0, beta, 1)] = x[std::slice(beta*(k+who*gamma), beta, 1)];
		bit_reverse_copy(aux, aux2);
		__1D_FFT_recursive(aux2);
		while (1)
		{
			if (m.try_lock()) {
				y[std::slice(beta*(k+who*gamma), beta, 1)] = aux2[std::slice(0, beta, 1)];
				m.unlock();
				break;
			}
		}
	}
}






void PFFT2_iterative(std::valarray<std::complex<double>> &x, std::valarray<std::complex<double>> &y, unsigned long columns, int thread_number)
{
	unsigned long samples, rows;
	std::valarray<std::complex<double>> aux;
	std::thread threads[thread_number];
	std::mutex m;

	samples = x.size();
	if (samples == 0) {
		return;
	}
	aux.resize(samples);
	y.resize(samples);
	rows = samples/columns;

	for (int k = 0; k < thread_number; k++) {
		threads[k] = std::thread(recursive_column_wrapper, std::ref(x), std::ref(aux), rows, columns, columns/thread_number, k, std::ref(m));
	}	
	for (int k = 0; k < thread_number; k++) {
		threads[k].join();
	}
	

	for (int k = 0; k < thread_number; k++) {
		threads[k] = std::thread(recursive_row_wrapper, std::ref(aux), std::ref(y), rows, columns, rows/thread_number, k, std::ref(m));
	}	
	for (int k = 0; k < thread_number; k++) {
		threads[k].join();
	}
}
