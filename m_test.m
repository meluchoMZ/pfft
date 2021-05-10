%% Executes the same test than my implementation to see the difference
%% Digital information processing, Computer engineering
%% Author: Miguel Blanco Godón

clear; close all;

file = fopen('./testbench/test_0.dat');
x = fscanf(file, '%d');
x = reshape(x, 256, 256);
disp('FFT2 execution on 256x256 matrix: ');
tic; fft2(x); toc;
file = fopen('./testbench/test_1.dat');
x = fscanf(file, '%d');
x = reshape(x, 512, 512);
disp('FFT2 execution on 512x512 matrix: ');
tic; fft2(x); toc;
file = fopen('./testbench/test_2.dat');
x = fscanf(file, '%d');
x = reshape(x, 1024, 1024);
disp('FFT2 execution on 1024x1024 matrix: ');
tic; fft2(x); toc;
file = fopen('./testbench/test_3.dat');
x = fscanf(file, '%d');
x = reshape(x, 4096, 4096);
disp('FFT2 execution on 4096x4096 matrix: ');
tic; fft2(x); toc;
file = fopen('./testbench/test_4.dat');
x = fscanf(file, '%d');
x = reshape(x, 8192, 8192);
disp('FFT2 execution on 8192x8192 matrix: ');
tic; fft2(x); toc;
