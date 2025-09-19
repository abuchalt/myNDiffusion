%% 2D Cartesian 2-Group Neutron Diffusion
% ------------------------------------------------------------------------------
% Adam Buchalter
%
% This is a power-iteration eigenvalue solver for the 2-dimensional multigroup 
% neutron diffusion equation in a cartesian slab geometry employing a pseudo-
% finite-volume discretization scheme
% ------------------------------------------------------------------------------
clear all; close all; clc;

%% Import Data
% ------------------------------------------------------------------------------
myCWD = pwd;
UO2 = readtable(fullfile(myCWD,'data\\2Group_UO2.csv'));
H2O = readtable(fullfile(myCWD,'data\\2Group_H2O.csv'));
% Cell Array of Materials
M = {UO2, H2O};

%% Computational Parameters
% ------------------------------------------------------------------------------
G = 2; % Number of energy groups

% Define mesh size
fprintf('Mesh size') % separate print and input b/c vscode extension
i_max = input('');
j_max = i_max;

% Define variables for power-iteration
residual = 1.0E5; % init residual
epsilon = 1.0E-16; % drive residual down to this value before terminating

%% Intelligently Find Smallest Diffusion Length for Non-Dimensionalization
% ------------------------------------------------------------------------------
% The smallest diffusion length is important because if our cells are larger
% than the diffusion length for a particular group then we cannot reliably
% capture transport phenomena

L = 50; % [cm] Arbitrarily Large Initial Guess
for mat = 1:numel(M) % For each material
    for group = 1:G % For each group
        thisL = sqrt(M{mat}.D(group)/M{mat}.Sigma_a(group)); % Characteristic Diffusion Length [cm]
        if thisL < L % If smaller, use as new normalization length
            L = thisL;
        end
    end
end

%% Nondimensional Domain Prep
% ------------------------------------------------------------------------------

% Define physical domain
A = 10; % fuel slab dimension [L]
B = 20; % reflector thickness [L]
size = A + 2*B;

% Unitless Constants
% Abar = A/L;
% Bbar = B/L;
sizebar = size/L;
% LSig_a = L*Sigma_a;
% LnuSig_f = L*nuSigma_f;

% Calculate step sizes
Deltax = size/(i_max-1)/2;
Deltay = Deltax;
Deltaxbar = sizebar/(i_max-1)/2;
Deltaybar = Deltaxbar;