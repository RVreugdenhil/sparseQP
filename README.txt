## Principal Component Hierarchy for Sparse Quadratic Programs

This folder contains the source code for the experiments in the paper.

## Installation guide

Run the install.m script to add __'/src'__  folders to the working directory of MATLAB.

## Utilisation guide

The package comes with 2 main functions: BR and DP, listed in the __'/src/pca_approach', which are described in the main paper. 

BR and DP require the following software
- YALMIP 
- MOSEK 

The test.m file contains a simple example of BR and DP.

The folder compute_eigenvalues contains the code to evaluate the eigenvalues of the dataset listed in the paper

The folder compare_sparse_regression contains the code for the comparison of BR and DP to the warm start and screening,

- The warm start and screening implementation are listed in the __'/src/warmstart' and __'/src/screening'.
- The implementation of the method by Beck and Eldar is listed in __'/src/beck'
- The method by Yuan et al. is listed in __'/src/KDD'
