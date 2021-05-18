# Principal Component Hierarchy for Sparse Quadratic Programs
Authors: Robbie Vreugdenhil, Viet Anh Nguyen, Arman Eftekhari, Peyman Mohajerin Esfahani 
Featured at ICML 2021

This folder contains the source code for the experiments in the paper.

## Contents
- [Installation guide](#Installation-guide)
- [Introduction](#introduction)
- [Utilisation guide](#Utilisation-guide)
- [DORO](#doro)
  - [CelebA](#celeba)
  - [CivilComments-Wilds](#civilcomments-wilds)


## Installation guide

Run the install.m script to add __'/src'__  folders to the working directory of MATLAB.
To run experiments on the datasets used in the original paper unzip the __'/data.zip'__ file in the repository

## Introduction
We propose two scalable optimization algorithms for sparse quadratic programming problems, coined as the "best response" and the "dual program", that can efficiently screen the potential indices of the nonzero elements of the original program. We show that the proposed methods are competitive with the existing screening methods in the current sparse regression literature, and it is particularly fast on instances with high number of measurements in experiments with both synthetic and real datasets.

## Utilisation guide

The package comes with 2 main functions: BR and DP, listed in the __'/src/pca_approach', which are described in Section 3 the main paper. 

BR and DP require the following software
- YALMIP 
- MOSEK 

The test.m file contains a simple example of BR and DP.

The folder __'/compute_eigenvalues' contains the code to evaluate the eigenvalues of the dataset listed in the paper

The folder __'/compare_sparse_regression' contains the code for the comparison of BR and DP to the warm start, screening, Beck alg 7 and the KDD method.

- The warm start and screening implementation are listed in the __'/src/warmstart' and __'/src/screening'.
- The implementation of the method by Beck and Eldar is listed in __'/src/beck'
- The method by Yuan et al. is listed in __'/src/KDD'
