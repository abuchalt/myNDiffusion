%% Steady-State Heat Diffusion in a Westinghouse PWR Fuel Rod 
% ------------------------------------------------------------------------------
% Adam Buchalter
%
% Desc.
% ------------------------------------------------------------------------------
clear all; close all; clc;

%% Data
% ------------------------------------------------------------------------------

% General
% Fuel Type: UO2
% Coolant + Moderator: H2O
% Structural Material: Zircaloy
P = 3411*1E6; % Thermal Output [MW_th -> W_th]

% Core Data
L = 366; % Active Height [cm]
D_eq = 337; % Equivalent Diameter [cm]
HDRatio = 1.09; % H/D
V = 32800*1E3; % Core Volume [L->cc]
M = 90200; % Fuel Mass [kg]

% Assembly Data
N_Rods = 50952; % Total Number of Fuel Rods

% Fuel Rod Data
pitch = 1.25; % Fuel Rod Pitch [cm]
myOD = 0.94; % Fuel Rod Outer Diameter [cm]
c = 0.0572; % Clad Thickness [cm]
s = 0.819; % Fuel-pellet Diameter [cm]
delta = 0.0082; % Pellet-clad gap [cm]

% Thermal Hydraulic Data
p = 155*1E5; % System Pressure [bar->Pa]
mdot = 62E6/3600; % Coolant Mass Flow Rate [kg/h->kg/s]
qprime = 178; % Average Linear Power Density [W/cm]
qprimeMax = 426; % Peak Linear Power Density [W/cm]
q2prime = 68.5; % Average Heat Flux [W/cm^2]
q2primeMax = 183; % Peak Heat Flux [W/cm^2]
T_in = 300+273.15; % Inlet Temperature [C->K]
T_out = 332+273.15; % Outlet Temperature [C->K]
T_c = 1788+273.15; % Centerline Fuel Temperature [C->K]