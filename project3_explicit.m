%% 2D Cartesian Heat Diffusion (EXPLICIT)
% ------------------------------------------------------------------------------
% Adam Buchalter
%
% This is an explicit solver for 2-dimensional Fourier-Biot heat diffusion in a 
% cartesian slab geometry using a cell-centered finite-volume discretization 
% scheme and Newton (convective) boundary conditions.
% ------------------------------------------------------------------------------
clear all; close all; clc;

%% Import Data
% ------------------------------------------------------------------------------
myCWD = pwd;
UO2 = readtable(fullfile(myCWD,'thermalData\\UO2.csv'));

% Physical params
totPwr = 3000; % Total Core Power [MW_{th}]
% fuelLength = 400; % Total Average Fuel Rod Length [cm]
% totLinPwr = 1E6*totPwr/fuelLength; % Total Linear Heat Generation [W/cm]

%% Computational Parameters
% ------------------------------------------------------------------------------
% Bulk Convective Fluid Temperature
T_infty = 100; % [K]

% Define mesh size
% fprintf('Quarter-mesh size') % separate print and input b/c vscode extension
% i_max = input('');
i_max = 65;
j_max = i_max;

%% Nondimensional Domain Prep
% ------------------------------------------------------------------------------

% Define physical domain
size = 125; % Domain Size [cm]
maxNodes = i_max*2 - 1; % Total number of nodes in domain
% fuelDim = 58; % Fuel Dimensions [Deltax] or [number of nodes]
% modDim = ceil((maxNodes-fuelDim)/2);

% Unitless Constants
T_r = T_infty; % [K]

% Calculate step sizes
Deltax = size/(i_max*2 - 1);
Deltay = Deltax;
% Deltaxbar = sizebar/(i_max*2 - 1);
% Deltaybar = Deltaxbar;

% CHANGE TO VERIFY FOURIER NUMBER
% if Deltax > L
%     fprintf('Warning: Grid Resolution is Larger than Diffusion Length')
% end

% Define x and y values in spatial domain (fully dimensional)
for i = 1:i_max
    for j = 1:j_max
        x(i,j) = Deltax*(i-1);
        y(i,j) = Deltay*(j-1);
    end
end

% Define x and y values in complete spatial domain
for i = 1:(2*i_max)-1
    for j = 1:(2*j_max)-1
        fullx(i,j) = Deltax*(i-1);
        fully(i,j) = Deltay*(j-1);
    end
end

% File Info
subfolder='results\\project2\\'+string((2*i_max)-1)+'x'+string((2*i_max)-1);
% subfolder='results\\project2\\'+string((2*i_max)-1)+'x'+string((2*i_max)-1)+'_4G';
mkdir(fullfile(myCWD,subfolder));

%% Build Coefficient Matrices
% ------------------------------------------------------------------------------

%% Solve
% ------------------------------------------------------------------------------

%% Visualize Results
% ------------------------------------------------------------------------------

%% Store Results
% ------------------------------------------------------------------------------

%% Functions
% ------------------------------------------------------------------------------

% Pointer mapping function for 2D->1D transform of solution
function k = pmap(i, j, i_max)
    k = i + (j-1)*i_max;
end

function symk = sympmap(i, j, i_max, j_max)
    symk = (j_max-j+1) + (i_max-i)*j_max;
end