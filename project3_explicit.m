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
H2O = readtable(fullfile(myCWD,'thermalData\\H2O.csv'));

% Cell Array of Materials
M = {UO2, H2O};

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

% Define time stepping
Deltat = 0.1; % [s]

%% Nondimensional Domain Prep
% ------------------------------------------------------------------------------

% Define physical domain
size = 125; % Domain Size [cm]
maxNodes = i_max*2 - 1; % Total number of nodes in domain
fuelDim = maxNodes; % Fuel Dimensions [Deltax] or [number of nodes]
modDim = ceil((maxNodes-fuelDim)/2);

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

% Specify materials in domain
quarterFuelDim = (fuelDim+1)/2;
for i = 1:i_max
    for j = 1:j_max
        k = pmap(i, j, i_max);
        if i <= quarterFuelDim && j_max-j+1 <= quarterFuelDim % because of reflection mapping j is flipped
            mat(k) = 1; % fuel
        else
            mat(k) = 2; % moderator
        end
    end
end

% Specify Internal Volumetric Heat Generation [W/m^3]
q3prime = zeros(i_max*j_max,1); % Init

% File Info
subfolder='results\\project2\\'+string((2*i_max)-1)+'x'+string((2*i_max)-1);
% subfolder='results\\project2\\'+string((2*i_max)-1)+'x'+string((2*i_max)-1)+'_4G';
mkdir(fullfile(myCWD,subfolder));

%% Build Coefficient Matrices
% ------------------------------------------------------------------------------
% Init coeff matrices
A = spalloc(i_max*j_max, i_max*j_max, 5*i_max*j_max); % Sparsely allocate Diffusion Operator with 5 bands
B = zeros(i_max*j_max, 1); % Allocate Heat Generation Vector

% Init Solution Variables (1D because we use pointer mapping)
T = ones(i_max*j_max,1);
% "Old" Solution for computing residual
T_old = ones(i_max*j_max,1);

% Define Coefficient Matrix
for i = 2:i_max-1
    for j = 2:j_max-1
        k = pmap(i, j, i_max);
        k_e = k + 1;
        k_w = k - 1;
        k_n = k + i_max;
        k_s = k - i_max;

        % matnum = mat(k); % what material
        thisk_k = M{mat(k)}.k;
        thisk_ke = M{mat(k_e)}.k;
        thisk_kw = M{mat(k_w)}.k;
        thisk_kn = M{mat(k_n)}.k;
        thisk_ks = M{mat(k_s)}.k;
        thisrhoc_p = M{mat(k)}.rhoc_p;

        nonDimFactorx = Deltat/(2.0*thisrhoc_p*Deltax^2);
        nonDimFactory = Deltat/(2.0*thisrhoc_p*Deltay^2);
        % pointer mapping goes row-by-row to assemble Coeff. Matrix
        A(k,k) = 1.0 - (nonDimFactorx*(thisk_ke+2.0*thisk_k+thisk_kw)) - (nonDimFactory*(thisk_kn+2.0*thisk_k+thisk_ks));
        A(k,k_e) = nonDimFactorx*(thisk_ke+thisk_k);
        A(k,k_w) = nonDimFactorx*(thisk_kw+thisk_k);
        A(k,k_n) = nonDimFactory*(thisk_kn+thisk_k);
        A(k,k_s) = nonDimFactory*(thisk_ks+thisk_k);

        B(k) = (Deltat/(thisrhoc_p*T_r))*q3prime(k);
    end
end
A = sparse(A); % Enforce Coeffs sparse

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