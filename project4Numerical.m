%% 1D Cylindrical Heat Diffusion
% ------------------------------------------------------------------------------
% Adam Buchalter
%
% words
% ------------------------------------------------------------------------------
clear all; close all; clc;

%% Import Data
% ------------------------------------------------------------------------------
myCWD = pwd;
fuel = readtable(fullfile(myCWD,'results\\project4\\thermalData\\fuel.csv'));
clad = readtable(fullfile(myCWD,'results\\project4\\thermalData\\clad.csv'));
gap = readtable(fullfile(myCWD,'results\\project4\\thermalData\\gap.csv'));
h_coolant = 3.5; % Coolant Convective Heat Transfer [W/cm^2.K]
q3prime_fuel = 347.20160; % Volumetric Heat Rate in Fuel [W/cm^3]
q_s = 66945.36034; % Total Heat Produced by a Fuel Element [W]

c = 0.0572; % Clad Thickness [cm]
s = 0.819; % Fuel-pellet Diameter [cm]
g = 0.0082; % Pellet-clad gap [cm]
R = s/2; % Fuel-pellet Radius [cm]

% Cell Array of Materials
M = {fuel, gap, clad};

%% Computational Parameters
% ------------------------------------------------------------------------------
% Bulk Convective Fluid Temperature
T_infty = 589.15; % [K]

% Define mesh size
nodesPerReg = 129;
N_fuel = nodesPerReg;
deltar_fuel = R/N_fuel;
N_gap = nodesPerReg;
deltar_gap = g/N_gap;
N_clad = nodesPerReg;
deltar_clad = c/N_clad;

i_max = N_fuel+N_gap+N_clad;

%% Domain Prep
% ------------------------------------------------------------------------------

% Specify materials in domain
for i = 1:i_max
    if i <= N_fuel
        mat(i) = 1;
        deltar(i) = deltar_fuel;
        r(i) = i*deltar_fuel - deltar_fuel/2;
    elseif i <= N_fuel+N_gap
        mat(i) = 2;
        deltar(i) = deltar_gap;
        r(i) = N_fuel*deltar_fuel + (i - N_fuel)*deltar_gap - deltar_gap/2;
    else
        mat(i) = 3;
        deltar(i) = deltar_clad;
        r(i) = N_fuel*deltar_fuel + N_gap*deltar_gap + (i - N_fuel - N_gap)*deltar_clad - deltar_clad/2;
    end
end

%% Build Coefficient Matrices
% ------------------------------------------------------------------------------
% Init coeff matrix
A = spalloc(i_max, i_max, 3*i_max); % Sparsely allocate Diffusion Operator with 3 bands
Q = zeros(i_max, 1); % Allocate Heat Generation Vector

% Define Coefficient Matrix
for i = 2:i_max-1
    thisk = M{mat(i)}.k;
    thisk_inner = M{mat(i-1)}.k;
    thisk_outer = M{mat(i+1)}.k;

    innerCoeff = 2.0 * ( r(i) - deltar(i)/2 ) * ( thisk_inner*deltar(i-1) + thisk*deltar(i) ) / ( r(i) * deltar(i) * ( deltar(i-1) + deltar(i) )^2 );
    outerCoeff = 2.0 * ( r(i) + deltar(i)/2 ) * ( thisk_outer*deltar(i+1) + thisk*deltar(i) ) / ( r(i) * deltar(i) * ( deltar(i+1) + deltar(i) )^2 );
    
    A(i,i-1) = -innerCoeff;
    A(i,i) = innerCoeff + outerCoeff;
    A(i,i+1) = - outerCoeff;
    
    if i <= N_fuel
        Q(i) = q3prime_fuel;
    else
        Q(i) = 0.0;
    end
end

% Innermost Boundary -- no-gradient
i = 1;

thisk = M{mat(i)}.k;
thisk_outer = M{mat(i+1)}.k;

outerCoeff = 2.0 * ( r(i) + deltar(i)/2 ) * ( thisk_outer*deltar(i+1) + thisk*deltar(i) ) / ( r(i) * deltar(i) * ( deltar(i+1) + deltar(i) )^2 );

A(i,i) = outerCoeff;
A(i,i+1) = -outerCoeff;

Q(i) = q3prime_fuel;

% Outermost Boundary -- Newton
i = i_max;

thisk = M{mat(i)}.k;
thisk_inner = M{mat(i-1)}.k;

innerCoeff = 2.0 * ( r(i) - deltar(i)/2 ) * ( thisk_inner*deltar(i-1) + thisk*deltar(i) ) / ( r(i) * deltar(i) * ( deltar(i-1) + deltar(i) )^2 );

A(i,i-1) = innerCoeff;
A(i,i) = -innerCoeff - h_coolant / ( r(i) * deltar(i) );

Q(i) = - h_coolant * T_infty / ( r(i) * deltar(i) );

A = sparse(A);

T = A\Q;

%% Visualize Results
% ------------------------------------------------------------------------------
figure(1);
plot(r, T);
ylabel('Temperature [K]');
xlabel('r [cm]');
title('Steady-State Solution for a Westinghouse PWR Fuel Pin');

%% Store Results
% ------------------------------------------------------------------------------
subsubfolder = [num2str(nodesPerReg),'nodes'];
plotOut = fullfile(myCWD,'results','project4',subsubfolder);
mkdir(plotOut);

% Save Figure
saveas(figure(1),fullfile(plotOut,'tempContour.jpg'));

% And Steady State Solution Matrix
save(fullfile(plotOut,'T.mat'), 'T');
save(fullfile(plotOut,'r.mat'), 'r');

%% Try Different h_coolant
% ------------------------------------------------------------------------------

h_coolant = 0.0192332; % Coolant Convective Heat Transfer [W/cm^2.K]
Tclim = 1850+273.15; % Clad Thermal Limit [C->K]
Tflim = 2800+273.15; % Fuel Thermal Limit [C->K]

% Outermost Boundary -- Newton
i = i_max;

thisk = M{mat(i)}.k;
thisk_inner = M{mat(i-1)}.k;

innerCoeff = 2.0 * ( r(i) - deltar(i)/2 ) * ( thisk_inner*deltar(i-1) + thisk*deltar(i) ) / ( r(i) * deltar(i) * ( deltar(i-1) + deltar(i) )^2 );

A(i,i-1) = innerCoeff;
A(i,i) = -innerCoeff - h_coolant / ( r(i) * deltar(i) );

Q(i) = - h_coolant * T_infty / ( r(i) * deltar(i) );

A = sparse(A);

T = A\Q;

Tfmax = T(1)
Tcmax = T(N_fuel+N_gap+1)

if Tfmax > Tflim
    fprintf('Fuel Limit Exceeded!')
end
if Tcmax > Tclim
    fprintf('Clad Limit Exceeded!')
end

%% Visualize Results
% ------------------------------------------------------------------------------
figure(2);
plot(r, T);
ylabel('Temperature [K]');
xlabel('r [cm]');
title('Meltdown Condition Profile');

%% Store Results
% ------------------------------------------------------------------------------
plotOut = fullfile(myCWD,'results','project4');
mkdir(plotOut);

% Save Figure
saveas(figure(2),fullfile(plotOut,'meltdownProfile.jpg'));