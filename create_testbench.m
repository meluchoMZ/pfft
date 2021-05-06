%% Creates several test files for checking the FFT computation
%% Digital information processing, Computer engineering
%% Author: Miguel Blanco Godón

clear; close all;

im = imread('./media/image_0.jpeg');
writematrix(im, './testbench/test_0.dat');
