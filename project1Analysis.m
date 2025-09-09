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
    phiPages(1:minN,1:minN,i) = phi(1:r:N,1:r:N);
    kList(i) = keff(1);
end

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
    midline = (N+1)/2;

    % Define x and y values in complete spatial domain
    Deltax = A/(N-1);
    fullx = zeros(N,1);
    for j = 1:N
        fullx(j) = Deltax*(j-1);
    end

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

% Init
errList = zeros(meshes, 1);

% Gather and Plot Errors in Solution Matrices
figure(3);
for i = 1:meshes
    N = meshSizes(i);
    subfolder=string(N)+'x'+string(N);
    myphi = matfile(fullfile(mydir,subfolder,'phiPlot.mat'));
    phi = reshape(myphi.phiref2Plot, N, N);
    midline = (N+1)/2;

    % Define x and y values in complete spatial domain
    Deltax = A/(N-1);
    fullx = zeros(N,1);
    for j = 1:N
        fullx(j) = Deltax*(j-1);
    end

    myErrProfile = phi(midline,:)-transpose(analytical(fullx, B/2, A, B));
    
    % errList(i) = sum(myErrProfile);

    plot(fullx, myErrProfile, 'DisplayName', subfolder);
    hold on;
end
hold off;
legend;
ylabel('Difference in Normalized Flux','interpreter','latex');
xlabel('x-Axis Distance [cm]','interpreter','latex');
title('Error in Centerline Flux Profile for Different Grid Spacings','interpreter','latex');
xlim([0 A])
filename = 'globalErrorFlux.jpg';
saveas(figure(3),fullfile(mydir,filename));

%% Calculate Global Errors
% ------------------------------------------------------------------------------

for i = 1:meshes
    N = meshSizes(i);
    subfolder=string(N)+'x'+string(N);
    myphi = matfile(fullfile(mydir,subfolder,'phiPlot.mat'));
    phi = reshape(myphi.phiref2Plot, N, N);
    midline = (N+1)/2;

    analyticalSurface = zeros(N,N);

    % Define x and y values in complete spatial domain
    Deltax = A/(N-1);
    Deltay = B/(N-1);
    fullx = zeros(N,1);
    for j = 1:N
        for k = 1:N
            x = Deltax*(j-1);
            y = Deltay*(k-1);
            analyticalSurface(j,k) = analytical(x, y, A, B);
        end
    end

    myErrSurface = phi-analyticalSurface;

    errList(i) = sum(myErrSurface,"all")/(N^2);
end

%% Graphically Estimate Order of Convergence by Global Error
% ------------------------------------------------------------------------------

r = 2;
p2 = log(norm(errList(1,1)-errList(2,1))/norm(errList(2,1)-errList(3,1)))/log(r);
fprintf('Order of Grid Convergence by Error: %g\n', p2);

% Demonstrate Order of Convergence Graphically
h = 1./meshSizes;
% Transform Into Expected Linear Behavior Space
h2 = h.^2;
% Line of Fit
coefficients = polyfit(h2, errList, 1);
xFit = linspace(min(h2), max(h2), 3);
yFit = polyval(coefficients , xFit);
% Inspect Closeness
figure(4);
plot(h2, errList, '.', 'MarkerSize', 15); % Plot Real Data
hold on;
plot(xFit, yFit, 'k--'); % Plot Fit
hold on;
yl.LabelHorizontalAlignment = 'center';
yl.Color = [.90 0 0];
hold off;
ylabel('Global Error in Flux Field','interpreter','latex');
xlabel('Square of Relative Grid Spacing $h^2$','interpreter','latex');
title('Visualize Convergence by Global Relative Error','interpreter','latex');
set ( gca, 'XDir', 'reverse' )
filename = 'GridConvergence_GlobalError.jpg';
saveas(figure(4),fullfile(mydir,filename));

%% Functions
% ------------------------------------------------------------------------------

% Analytical Flux Function
function phi = analytical(x, y, A, B)
    phi = sin(pi.*x/A)*sin(pi.*y/A)/(4*A*B/pi^2);
end