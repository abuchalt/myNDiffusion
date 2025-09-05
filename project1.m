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
%          o  <---J_w---   o   ---J_e--->  o        \vec{J} = -D \grad\phi
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
A = 1.0;
B = 1.0;

% Input parameters
Sigma_tr = 3.62E-2; % Macroscopic Transport Cross-Section [cm^-1]
Sigma_a = 0.1532; % Macroscopic Absorption Cross-Section [cm^-1]
nuSigma_f = 0.1570; % Product of Macroscopic Fission Cross-Section and Neutrons per Fission [cm^-1]
% Calculated Constant Parameters
D = 1.0/(3.0*(Sigma_a+Sigma_tr)); % Diffusion Coefficient [cm]
L = sqrt(D/Sigma_a); % Characteristic Diffusion Length [cm]
% Unitless Constants
Abar = A/L;
Bbar = B/L;
LSig_a = L*Sigma_a;
LnuSig_f = L*nuSigma_f;

% Define mesh size
fprintf('Maximum number of points in x-direction') % separate print and input
                                                   % b/c vscode extension
i_max = input('');
fprintf('Maximum number of points in y-direction')
j_max = input('');

% Calculate step sizes
Deltax = A/(i_max-1);
Deltay = B/(j_max-1);
Deltaxbar = Abar/(i_max-1);
Deltaybar = Bbar/(j_max-1);

% Define variables for power-iteration
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
myCWD = pwd;
subfolder='results\\project1\\'+string(i_max)+'x'+string(j_max);
mkdir(fullfile(myCWD,subfolder));

%% Script
% ------------------------------------------------------------------------------
% Init coeff matrices
M = spalloc(i_max*j_max, i_max*j_max, 5*i_max*j_max); % Sparsely allocate transport matrix with 5 bands
F = spalloc(i_max*j_max, i_max*j_max, 1*i_max*j_max); % Sparsely allocate fission matrix with 1 band

% Init Solution Variables (1D because we use pointer mapping)
psi = ones(i_max*j_max,1);
k = 1;
% "Old" Solution for computing residual
psi_old = ones(i_max*j_max,1);
k_old = 1;

% Define Coefficient Matrix
for i = 2:i_max-1
    for j = 2:j_max-1
        k = pmap(i, j, i_max);
        k_e = k + 1;
        k_w = k - 1;
        k_n = k + i_max;
        k_s = k - i_max;

        % pointer mapping goes row-by-row to assemble Coeff. Matrix
        M(k,k) = LSig_a*(1.0 + (2.0/Deltaxbar^2) + (2.0/Deltaybar^2));
        M(k,k_e) = (-1.0*LSig_a/Deltaxbar^2);
        M(k,k_w) = (-1.0*LSig_a/Deltaxbar^2);
        M(k,k_n) = (-1.0*LSig_a/Deltaybar^2);
        M(k,k_s) = (-1.0*LSig_a/Deltaybar^2);

        F(k,k) = (1.0*LnuSig_f/k);
    end
end
M = sparse(M); % Enforce Coeffs sparse
F = sparse(F);

%% Functions
% ------------------------------------------------------------------------------

% Pointer mapping function for 2D->1D transform of solution
function k = pmap(i, j, i_max)
    k = i + (j-1)*i_max;
end