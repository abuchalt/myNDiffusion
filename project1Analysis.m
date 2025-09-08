%% Project 1 Analysis
% ------------------------------------------------------------------------------
% Adam Buchalter
%
% This is a script for:
% - Estimating and visualizing convergence of truncation error for the 2-D
%   cartesian slab neutron transport problem
% - Visualizing flux line profiles
% - Comparing numerical solution to analytical solution
% ------------------------------------------------------------------------------
clear all; close all; clc;

%% User Input
% ------------------------------------------------------------------------------

% Params
meshSizes = [257, 129, 65, 33, 17, 9, 5];

% File Info
myCWD = pwd;
% subfolder='results\\project1';
subfolder='results\\project1Sym';
mydir=fullfile(myCWD,subfolder);

%% Import Data
% ------------------------------------------------------------------------------

% Quick Maths
meshSizes = sort(meshSizes);
minN = meshSizes(1);
meshes = size(meshSizes, 2);
assert(meshes>=3, 'Need at Least 3 Grids')
midline = (minN+1)/2;

% Init
phiPages = zeros(minN,minN,meshes);
kList = zeros(meshes, 1);

% Evaluate Error only at coincident points of uniform meshes
% function E = getError(f_fine, f_coarse)
%     Nfine = (size(f_fine,1)-1)
%     Ncoarse = (size(f_coarse,1)-1)
%     r = round(Nfine/Ncoarse);
%     f_fineresize = f_fine(1:r:Nfine);
%     E = size(f_fineresize,1);
% end

% Store Resized Solution Matrices
for i = 1:meshes
    N = meshSizes(i);
    r = (N-1)/(minN-1);
    assert(floor(r)==r, 'Grid Refinement Ratio is Non-Integer');
    subfolder=string(N)+'x'+string(N);
    myphi = matfile(fullfile(mydir,subfolder,'phiPlot.mat'));
    myk = matfile(fullfile(mydir,subfolder,'keff.mat'));
    phi = reshape(myphi.phiref2Plot, N, N);
    keff = myk.keff;
    phiPages(1:minN,1:minN,i) = phi(1:r:N,1:r:N)/norm(phi(1:r:N,1:r:N));
    kList(i) = keff(1);
end

%% Estimate Order of Convergence by Flux at Center
% ------------------------------------------------------------------------------

% Init
f = zeros(meshes, 1);

% Evaluate Order of Convergence Directly
for i = 1:meshes
    f(i, 1) = phiPages(midline,midline,i);
end
r = 2;
p2 = log(norm(f(3,1)-f(2,1))/norm(f(2,1)-f(1,1)))/log(r);
fprintf('Order of Grid Convergence by Flux at Center: %g\n', p2);

%% Graphically Estimate Order of Convergence by keff
% ------------------------------------------------------------------------------

% Init
f = zeros(meshes, 1);

% Evaluate Order of Convergence Directly
for i = 1:meshes
    f(i, 1) = kList(i);
end
r = 2;
p2 = log(norm(f(1,1)-f(2,1))/norm(f(2,1)-f(3,1)))/log(r);
fprintf('Order of Grid Convergence by k_eff: %g\n', p2);

% Demonstrate Order of Convergence Graphically
h = 1./meshSizes;
% Transform Into Expected Linear Behavior Space
h2 = h.^2;
% Line of Fit
coefficients = polyfit(h2, kList, 1);
xFit = linspace(min(h2), max(h2), 3);
yFit = polyval(coefficients , xFit);
% Inspect Closeness
figure(1);
plot(h2, kList, '.', 'MarkerSize', 15); % Plot Real Data
hold on;
plot(xFit, yFit, 'k--'); % Plot Fit
hold on;
yl = yline(1,'--','Analytical $k_{eff}$','interpreter','latex');
yl.LabelHorizontalAlignment = 'center';
yl.Color = [.90 0 0];
hold off;
ylabel('Predicted Neutron Multiplication Factor $k_{eff}$','interpreter','latex');
xlabel('Square of Relative Grid Spacing $h^2$','interpreter','latex');
title('Estimation of Order of Convergence by $k_{eff}$','interpreter','latex');
set ( gca, 'XDir', 'reverse' )
filename = 'GridConvergence_keff.jpg';
saveas(figure(1),fullfile(mydir,filename));

%% Plot Centerline Profiles
% ------------------------------------------------------------------------------

% Input parameters
Sigma_tr = 3.62E-2; % Macroscopic Transport Cross-Section [cm^-1]
Sigma_a = 0.1532; % Macroscopic Absorption Cross-Section [cm^-1]
nuSigma_f = 0.1570; % Product of Macroscopic Fission Cross-Section and Neutrons per Fission [cm^-1]
% Calculated Constant Parameters
D = 1.0/(3.0*(Sigma_a+Sigma_tr)); % Diffusion Coefficient [cm]
L = sqrt(D/Sigma_a); % Characteristic Diffusion Length [cm]

% Analytical Critical Dimension [cm]
critDim = pi * sqrt(2.0*D/(nuSigma_f-Sigma_a));

% Define physical domain
A = critDim; % x-Length [cm]
B = critDim; % y-Width [cm]

% Gather and Plot Solution Matrices
figure(2);
for i = 1:meshes
    N = meshSizes(i);
    subfolder=string(N)+'x'+string(N);
    myphi = matfile(fullfile(mydir,subfolder,'phiPlot.mat'));
    phi = reshape(myphi.phiref2Plot, N, N);

    % Define x and y values in complete spatial domain
    Deltax = A/(N-1);
    for j = 1:N
        fullx(j) = Deltax*(j-1);
    end
    phi = phi/(Deltax^2); % Renormalize by Mesh Spacing to target same graphical height i.e. multiply by cell areas per cm^2
    % output of program is (unitless power)/(cell area), and cell area changes with grid refinement... consider building this into diffusion program to get a result in W/cm^2

    plot(fullx, phi(midline,:), 'DisplayName', subfolder);
    hold on;
end
hold off;
legend;
ylabel('Normalized Neutron Flux $\phi$','interpreter','latex');
xlabel('x-Axis Distance [cm]','interpreter','latex');
title('Centerline Flux Profile for Different Grid Spacings','interpreter','latex');
xlim([0 A])
filename = 'lineProfileFlux.jpg';
saveas(figure(2),fullfile(mydir,filename));