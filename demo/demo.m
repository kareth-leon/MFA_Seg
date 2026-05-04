%==========================================================================
% Script Name : demo.m
%
% Description :
%   Implements a demo of the paper "Bayesian Multifractal Image
%   Segmentation" for the MRW with K=2 with multifractal values
%   theta1=-0.02 and theta2 = -0.08 (Figure 4-a)
%
% Dependencies :
%   - MATLAB R2022b or later
%   - Statistics and Machine Learning Toolbox
%   - Some functions of the MFA Toolbox found 
%     in https://www.irit.fr/~Herwig.Wendt/software.html#wlbmf are
%     incorporated in src_MGA_SEG
%
% Reference :
%   León-López, Kareth M., et al. "Bayesian Multifractal Image Segmentation."
%   IEEE Transactions on Image Processing 34 (2025): 8500-8510.
% 
% Authors :
%   Kareth LEON
%   Herwig WENDT
%   IRIT Laboratory, Université de Toulouse, CNRS, 
%   INP-Toulouse, UT3, UT2, TéSA, Toulouse, France.
%
% Date :
%   May 2026
%
% License :
%   MIT License
%==========================================================================
clear
close all 
clc
addpath('data\')
addpath('tools\')
addpath('trans_matrices\')

%--------------------------------------------------------------------------
% Load data
load('70_30_100_2_256x256_mrw2D_TWO_PARTS_rtiox2_10.mat')
vv0 = 100; % realization to be used
X = data(:,:,vv0);
%--------------------------------------------------------------------------
% Procedure
parameters_initialization
trans_mat % down/up matrices
%--------------------------------------------------------------------------
gibbs_sampler
%--------------------------------------------------------------------------
% Performance
performance_eval
%--------------------------------------------------------------------------
