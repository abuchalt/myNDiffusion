%% 2D Slab Neutron Diffusion
% ------------------------------------------------------------------------------
% Adam Buchalter
%
% This is a power-iteration eigenvalue solver for the 2-dimensional neutron 
% diffusion equation in a cartesian slab geometry employing a pseudo-finite-
% volume discretization scheme
% ------------------------------------------------------------------------------
clear all; close all; clc;

% ==============================================================================
%   Our governing equation:
%
%       -D \nabla^2 \phi + \Sigma_a \phi = \frac{1}{k} \nu \Sigma_f \phi
%
% ==============================================================================
%   Behold, an internal node:
%
%          o       |       o       |       o
%                  |       ^       |        
%                  |       |       |        
%          --------|------J_n------|--------
%                  |       |       |        
%                  |               |        
%          o  <---J_w---   o   ---J_e--->  o        \vec{J} = D \grad\phi
%                  |               |        
%                  |       |       |        
%          --------|------J_s------|--------
%                  |       |       |        
%                  |       V       |        
%          o       |       o       |       o
%
% ==============================================================================

%% Computational Parameters
% ------------------------------------------------------------------------------

% Define physical domain
h = 1.0;
w = 1.0;

% Define mesh size
fprintf('Maximum number of points in x-direction') % separate print and input
                                                   % b/c vscode extension
i_max = input('');
fprintf('Maximum number of points in y-direction')
j_max = input('');

% Calculate step sizes
Deltax = w/(i_max-1);
Deltay = h/(j_max-1);

% Input parameters

% And define variables for power-iteration
residual = 1.0E5; % init residual
epsilon = 1.0E-12; % drive residual down to this value before terminating

% Define x and y values in spatial domain
for i = 1:i_max
    for j = 1:j_max
        x(i,j) = Deltax*(i-1);
        y(i,j) = Deltay*(j-1);
    end
end

% File Info
myCWD = pwd
subfolder='results\\project1\\'+string(i_max)+'x'+string(j_max);
mkdir(fullfile(myCWD,subfolder));

%% Script
% ------------------------------------------------------------------------------

%% Functions
% ------------------------------------------------------------------------------

function k = pmap(i, j, i_max)
    k = i + (j-1)*i_max;
end