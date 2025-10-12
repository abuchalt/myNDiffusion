%% Project 3 Analysis
% ------------------------------------------------------------------------------
% Adam Buchalter
%
% This is a script for:
% - Estimating, visualizing, and comparing convergence of truncation error
%   between explicit and implicit formulations of the 2-D cartesian slab
%   heat diffusion problem in steady-state and transient
% - Comparing computational efficiency (CPU-time) among the explicit and
%   implicit formulations
% - Visualizing temperature line profiles
% ------------------------------------------------------------------------------
clear all; close all; clc;

%% User Input
% ------------------------------------------------------------------------------

% Params
meshSizes = [65, 33, 17];
timeSteps = [0.4, 0.2, 0.1];
myCWD = pwd;

% Quick Maths
meshSizes = sort(meshSizes);
minN = meshSizes(1);
meshes = size(meshSizes, 2);
assert(meshes>=3, 'Need at Least 3 Grids')
midline = (minN+1)/2;

expFolder = 'results\\project3exp';
impFolder = 'results\\project3exp';

% ------------------------------------------------------------------------------