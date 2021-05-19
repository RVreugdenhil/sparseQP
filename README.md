# Principal Component Hierarchy for Sparse Quadratic Programs
Authors: Robbie Vreugdenhil, Viet Anh Nguyen, Arman Eftekhari, Peyman Mohajerin Esfahani 
Featured at ICML 2021

This folder contains the source code for the experiments in the paper.

## Contents
- [Installation guide](#Installation-guide)
- [Introduction](#introduction)
  - [Eigenvalues Real Data](#eigenvalues-real-data)
- [Best Response and Dual Program](#BR-and-DP)
  - [Example](#Example)
  - [Sparse Regression](#sparse-regression)


## Installation guide

Run the install.m script to add __'/src'__  folders to the working directory of MATLAB.
To run experiments on the datasets used in the original paper unzip the __'/data.zip'__ file in the repository

## Introduction
We propose two scalable optimization algorithms for sparse quadratic programming problems, coined as the "best response" and the "dual program", that can efficiently screen the potential indices of the nonzero elements of the original program. We show that the proposed methods are competitive with the existing screening methods in the current sparse regression literature, and it is particularly fast on instances with high number of measurements in experiments with both synthetic and real datasets.

### Eigenvalues Real Data

The folder __'/compute_eigenvalues'__ contains the code to evaluate the eigenvalues of the datasets listed in the paper. Running __'/compute_eigenvalues/analyse_eig_decay.m'__ will give the same results as listed in Section 1 in the main paper using the datasets in __'/data.zip'__ . These results can also be found in the  __'/compute_eigenvalues/results'__ folder.

## Best Response and Dual Program

In Section 4, we conduct numerical experiments on both synthetic and real datasets. Here we describe how to run the experiments

The package comes with 2 main functions: BR (Best Response) and DP (Dual Program), listed in the __'/src/pca_approach'__, which are described in Section 3 the main paper. 

BR and DP require the following software
- YALMIP 
- MOSEK 

### Example 

The test.m file contains a simple example of BR and DP on a synthetic dataset.

### Sparse Regression

The folder __'/compare_sparse_regression'__ contains the code for the comparison of BR and DP to the warm start, screening, Beck alg 7 and the KDD method. 

- The warm start and screening implementation are listed in the __'/src/warmstart'__ and __'/src/screening'__.
- The implementation of the method by Beck and Eldar is listed in __'/src/beck'__
- The method by Yuan et al. is listed in __'/src/KDD'__

The __'/compare_sparse_regression/different_k'__ folder gives a comparison between BR and DP for different values of k, where k is described in the paper. 

The __'/compare_sparse_regression/syn_constant_size'__ folder gives a comparison of all different methods for a synthetic dataset of constant size N and n but the other parameters are variable (e.g. \rho, s, \eta)

The __'/compare_sparse_regression/syn_different_size'__ folder gives a comparison of all different methods for a synthetic dataset with variable N and n but the other parameters are fixed (e.g. \rho, s, \eta)

The __'/compare_sparse_regression/real'__ folder gives a comparison of all different methods over 5 different real datasets. 

