close all; warning off;
clear;
clc;
addpath(genpath('LIB'));
addpath(genpath('ClusteringMeasure'));
load('./Data_Ting.mat')
num_cluster = max(unique(true_labs));
alpha = 1e-4;
lambda = 1e-3;
delta = 1e-4;
% run
[U] = DcGSE(in_X, true_labs, alpha, lambda, delta);