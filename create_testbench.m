%% Creates several test files for checking the FFT computation
%% Digital information processing, Computer engineering
%% Author: Miguel Blanco Godón

clear; close all;

im = imread('./media/image_0.jpeg');
% 256x256 submatrix from image
x1 = reshape(im(1:256, 1:256), [], 1);
% 512x512 submatrix from image
x2 = reshape(im(1:512, 1:512), [], 1);
% 1024x1024 submatrix from image
x3 = reshape(im(1:1024, 1:1024), [],1);
% 4096x4096 random generated matrix
x4 = floor(255.*rand(4096*4096,1));
% 8192x8192 random generated matrix
x5 = floor(255.*rand(8192*8192, 1));
% 2^14x2^14 random generated matrix
x6 = floor(255.*rand(2^14 * 2^14, 1));

writematrix(x1, './testbench/test_0.dat');
writematrix(x2, './testbench/test_1.dat');
writematrix(x3, './testbench/test_2.dat');
writematrix(x4, './testbench/test_3.dat');
writematrix(x5, './testbench/test_4.dat');
writematrix(x6, './testbench/test_5.dat');
