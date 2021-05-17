clc;clear all;close all;

cd 'solver'
mex -O -lmwlapack -lmwblas QuadFractionalProgrammingCOO.c;
cd ..

