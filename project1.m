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
epsilon = 1.0E-16; % drive residual down to this value before terminating

% Define x and y values in spatial domain (fully dimensional)
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

%% Build Coefficient Matrices
% ------------------------------------------------------------------------------
% Init coeff matrices
M = spalloc(i_max*j_max, i_max*j_max, 5*i_max*j_max); % Sparsely allocate transport matrix with 5 bands
F = spalloc(i_max*j_max, i_max*j_max, 1*i_max*j_max); % Sparsely allocate fission matrix with 1 band

% Init Solution Variables (1D because we use pointer mapping)
phi = ones(i_max*j_max,1);
keff = 1;
% "Old" Solution for computing residual
phi_old = ones(i_max*j_max,1);
keff_old = 1;

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

        F(k,k) = (1.0*LnuSig_f);
    end
end
M = sparse(M); % Enforce Coeffs sparse
F = sparse(F);

% Apply BCs on initial guess
% Left BC
for i = 1:1
    for j = 1:j_max
        k = pmap(i,j,i_max);
        % k_e = k + 1;
        % k_ee = k + 2;
        phi(k) = 0; % Zero-flux
        M(k,k) = 1;
        F(k,k) = 1;
    end
end
% Right BC
for i = i_max:i_max
    for j = 1:j_max
        k = pmap(i,j,i_max);
        % k_w = k - 1;
        % k_ww = k - 2;
        phi(k) = 0; % Zero-flux
        M(k,k) = 1;
        F(k,k) = 1;
    end
end
% Bottom BC
for i = 1:i_max
    for j = 1:1
        k = pmap(i,j,i_max);
        % k_n = k + i_max;
        % k_nn = k + 2*i_max;
        phi(k) = 0; % Zero-flux
        M(k,k) = 1;
        F(k,k) = 1;
    end
end
% Top BC
for i = 1:i_max
    for j = j_max:j_max
        k = pmap(i,j,i_max);
        % k_s = k - i_max;
        % k_ss = k - 2*i_max;
        phi(k) = 0; % Zero-flux
        M(k,k) = 1;
        F(k,k) = 1;
    end
end

%% Power-Iteration Script
% ------------------------------------------------------------------------------
fprintf('test')
% Compute evolution operator initially to minimize work in loop
% Amat = inv(M)*F % Slowww
Amat = M\F;

% Init iteration vars
tTot = 0;
iter = 0;

% Begin iteration
while (residual > epsilon)
% while (iter<2)

    tStart = tic;

    % Track previous source vector, flux vector, and k
    phi_old = phi;
    keff_old = keff;

    % Solve new guess of Phi
    phi = Amat * phi; % Evolve Flux
    phi = phi/norm(phi); % Normalize
    phiT = transpose(phi);
    keff = phiT * (Amat * phi); % Search Dominant Eigenvalue

    % Compute the new residual
    residual = norm(phi-phi_old);
    residual = residual/(i_max*j_max); % Normalize for DOF


    % Plot solution
    if mod(iter,100) == 0
        phiPlot = reshape(phi, i_max, j_max);
        figure(1);
        % Plot flux surface
        surf(x,y,phiPlot);
        ylabel('y');
        xlabel('x');
        title('Flux Surface');
        drawnow;
    end

    tTot = tTot + toc(tStart);

    fprintf(1,'iter = %i, residual = %g\n',iter,log10(residual));
    iter = iter + 1;
end

%% Visualize Results
% ------------------------------------------------------------------------------

phiPlot = reshape(phi, i_max, j_max);
figure(1);
% Plot flux surface
surf(x,y,phiPlot);
ylabel('y');
xlabel('x');
title('Flux Surface');

fprintf(1,'keff = %f\n',keff);

%% Store Results
% ------------------------------------------------------------------------------

saveas(figure(1),fullfile(myCWD,subfolder,'fluxContour.jpg'));

% And Solution Matrix
save(fullfile(myCWD,subfolder,'phi.mat'), 'phi');
save(fullfile(myCWD,subfolder,'keff.mat'), 'keff');

% And Timing Info
tAvg = tTot/iter;
fid = fopen(fullfile(myCWD,subfolder,'time.txt'),'wt');
fprintf(fid, 'Total CPU-time: %s s\nAverage Time per Iteration: %s s\n', string(tTot), string(tAvg));
fclose(fid);

%% Functions
% ------------------------------------------------------------------------------

% Pointer mapping function for 2D->1D transform of solution
function k = pmap(i, j, i_max)
    k = i + (j-1)*i_max;
end