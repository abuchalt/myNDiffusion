%% 2D Cartesian Heat Diffusion (IMPLICIT)
% ------------------------------------------------------------------------------
% Adam Buchalter
%
% This is an implicit solver for 2-dimensional Fourier-Biot heat diffusion in a 
% cartesian slab geometry using a cell-centered finite-volume discretization 
% scheme and Newton (convective) boundary conditions.
% ------------------------------------------------------------------------------
clear all; close all; clc;

%% Import Data
% ------------------------------------------------------------------------------
myCWD = pwd;
% UO2 = readtable(fullfile(myCWD,'thermalData\\UO2.csv'));
% H2O = readtable(fullfile(myCWD,'thermalData\\H2O.csv'));
myMat = readtable(fullfile(myCWD,'thermalData\\example.csv')); % Example Thermal Properties

% Cell Array of Materials
% M = {UO2, H2O};
M = {myMat};

% Physical params
% totPwr = 3000; % Total Core Power [MW_{th}]
% fuelLength = 400; % Total Average Fuel Rod Length [cm]
% totLinPwr = 1E6*totPwr/fuelLength; % Total Linear Heat Generation [W/cm]

%% Computational Parameters
% ------------------------------------------------------------------------------
% Bulk Convective Fluid Temperature
T_infty = 20; % [K]

% Define mesh size
% fprintf('Quarter-mesh size') % separate print and input b/c vscode extension
% i_max = input('');
i_max = 17;
j_max = i_max;

% Define time stepping
Deltat = 0.1; % [s]

%% Nondimensional Domain Prep
% ------------------------------------------------------------------------------

% Define physical domain
size = 12.5; % Domain Size [cm]
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

% VERIFY FOURIER NUMBER
% Fo = 0; % Arbitrarily Small Initial Guess
% for thisMat = 1:numel(M) % For each material
%     thisk = M{thisMat}.k;
%     thisrhoc_p = M{thisMat}.rhoc_p;
%     thish = M{thisMat}.h;

%     thisFo = thisk*Deltat/(thisrhoc_p*Deltax^2);
    
%     if thisFo > Fo % If larger, use as new fourier number
%         Fo = thisFo;
%     end
% end
% fprintf(1,'Fo = %f\n',Fo);

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
for i = 1:i_max
    for j = 1:j_max
        k = pmap(i, j, i_max);
        mat(k) = 1; % Only One Material
    end
end

% quarterFuelDim = (fuelDim+1)/2;
% for i = 1:i_max
%     for j = 1:j_max
%         k = pmap(i, j, i_max);
%         if i <= quarterFuelDim && j_max-j+1 <= quarterFuelDim % because of reflection mapping j is flipped
%             mat(k) = 1; % fuel
%         else
%             mat(k) = 2; % moderator
%         end
%     end
% end

% Specify Internal Volumetric Heat Generation [W/cm^3]
q3prime = ones(i_max*j_max,1); % Uniform heating as in examples

% File Info
subfolder='results\\project3imp\\'+string((2*i_max)-1)+'x'+string((2*i_max)-1);
% subfolder='results\\project2\\'+string((2*i_max)-1)+'x'+string((2*i_max)-1)+'_4G';
mkdir(fullfile(myCWD,subfolder));

%% Build Coefficient Matrices
% ------------------------------------------------------------------------------
% Init coeff matrices
A = spalloc(i_max*j_max, i_max*j_max, 5*i_max*j_max); % Sparsely allocate Diffusion Operator with 5 bands
Q = zeros(i_max*j_max, 1); % Allocate Heat Source Vector

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
        %
        nonDimFactorx = Deltat/(2.0*thisrhoc_p*Deltax^2);
        nonDimFactory = Deltat/(2.0*thisrhoc_p*Deltay^2);
        % pointer mapping goes row-by-row to assemble Coeff. Matrix
        A(k,k) = 1.0 + (nonDimFactorx*(thisk_ke+2.0*thisk_k+thisk_kw)) + (nonDimFactory*(thisk_kn+2.0*thisk_k+thisk_ks));
        A(k,k_e) = -nonDimFactorx*(thisk_ke+thisk_k);
        A(k,k_w) = -nonDimFactorx*(thisk_kw+thisk_k);
        A(k,k_n) = -nonDimFactory*(thisk_kn+thisk_k);
        A(k,k_s) = -nonDimFactory*(thisk_ks+thisk_k);
        %
        Q(k) = (Deltat/(thisrhoc_p*T_r))*q3prime(k);
    end
end
A = sparse(A); % Enforce Coeffs sparse

% Apply BCs
% Left BC
for i = 1:1
    for j = 2:j_max-1 % Avoid corners
        k = pmap(i,j,i_max);
        k_e = k + 1;
        k_n = k + i_max;
        k_s = k - i_max;
        k_w = sympmap(i,j,i_max,j_max)-i_max;

        % Treat as normal internal node with special mapping
        % matnum = mat(k); % what material
        thisk_k = M{mat(k)}.k;
        thisk_ke = M{mat(k_e)}.k;
        thisk_kw = M{mat(k_w)}.k;
        thisk_kn = M{mat(k_n)}.k;
        thisk_ks = M{mat(k_s)}.k;
        thisrhoc_p = M{mat(k)}.rhoc_p;
        %
        nonDimFactorx = Deltat/(2.0*thisrhoc_p*Deltax^2);
        nonDimFactory = Deltat/(2.0*thisrhoc_p*Deltay^2);
        % pointer mapping goes row-by-row to assemble Coeff. Matrix
        A(k,k) = 1.0 + (nonDimFactorx*(thisk_ke+2.0*thisk_k+thisk_kw)) + (nonDimFactory*(thisk_kn+2.0*thisk_k+thisk_ks));
        A(k,k_e) = A(k,k_e) - nonDimFactorx*(thisk_ke+thisk_k);
        A(k,k_w) = A(k,k_w) - nonDimFactorx*(thisk_kw+thisk_k);
        A(k,k_n) = A(k,k_n) - nonDimFactory*(thisk_kn+thisk_k);
        A(k,k_s) = A(k,k_s) - nonDimFactory*(thisk_ks+thisk_k);
        %
        Q(k) = (Deltat/(thisrhoc_p*T_r))*q3prime(k);
    end
end
% Right BC -- Newton BC
for i = i_max:i_max
    for j = 2:j_max-1 % Avoid Corners
        k = pmap(i,j,i_max);
        k_w = k - 1;
        k_n = k + i_max;
        k_s = k - i_max;
        
        thisk_k = M{mat(k)}.k;
        thish_ke = M{1}.h; % In complete solver, pull in h from the east
        thisk_kw = M{mat(k_w)}.k;
        thisk_kn = M{mat(k_n)}.k;
        thisk_ks = M{mat(k_s)}.k;
        thisrhoc_p = M{mat(k)}.rhoc_p;

        nonDimFactorx = Deltat/(2.0*thisrhoc_p*Deltax^2);
        nonDimFactory = Deltat/(2.0*thisrhoc_p*Deltay^2);
        % pointer mapping goes row-by-row to assemble Coeff. Matrix
        A(k,k) = 1.0 + (nonDimFactorx*(2.0*thish_ke*Deltax+thisk_k+thisk_kw)) + (nonDimFactory*(thisk_kn+2.0*thisk_k+thisk_ks));
        A(k,k_w) = A(k,k_w) - nonDimFactorx*(thisk_kw+thisk_k);
        A(k,k_n) = A(k,k_n) - nonDimFactory*(thisk_kn+thisk_k);
        A(k,k_s) = A(k,k_s) - nonDimFactory*(thisk_ks+thisk_k);
        %
        Q(k) = (Deltat/(thisrhoc_p*T_r))*q3prime(k) + (thish_ke*Deltat*T_infty/(thisrhoc_p*Deltax*T_r));
    end
end
% Bottom BC -- Newton BC
for i = 2:i_max-1 % Avoid Corners
    for j = 1:1
        k = pmap(i,j,i_max);
        k_e = k + 1;
        k_w = k - 1;
        k_n = k + i_max;
        
        thisk_k = M{mat(k)}.k;
        thisk_ke = M{mat(k_e)}.k;
        thisk_kw = M{mat(k_w)}.k;
        thisk_kn = M{mat(k_n)}.k;
        thish_ks = M{1}.h; % In complete solver, pull in h from the south
        thisrhoc_p = M{mat(k)}.rhoc_p;

        nonDimFactorx = Deltat/(2.0*thisrhoc_p*Deltax^2);
        nonDimFactory = Deltat/(2.0*thisrhoc_p*Deltay^2);
        % pointer mapping goes row-by-row to assemble Coeff. Matrix
        A(k,k) = 1.0 + (nonDimFactorx*(thisk_ke+2.0*thisk_k+thisk_kw)) + (nonDimFactory*(thisk_kn+thisk_k+2.0*thish_ks*Deltay));
        A(k,k_e) = A(k,k_e) - nonDimFactorx*(thisk_ke+thisk_k);
        A(k,k_w) = A(k,k_w) - nonDimFactorx*(thisk_kw+thisk_k);
        A(k,k_n) = A(k,k_n) - nonDimFactory*(thisk_kn+thisk_k);
        %
        Q(k) = (Deltat/(thisrhoc_p*T_r))*q3prime(k) + (thish_ks*Deltat*T_infty/(thisrhoc_p*Deltay*T_r));
    end
end
% Top BC
for i = 2:i_max-1 % Avoid corners
    for j = j_max:j_max
        k = pmap(i,j,i_max);
        k_e = k + 1;
        k_w = k - 1;
        k_s = k - i_max;
        k_n = sympmap(i,j,i_max,j_max)+1;

        % Treat as normal internal node with special mapping
        % matnum = mat(k); % what material
        thisk_k = M{mat(k)}.k;
        thisk_ke = M{mat(k_e)}.k;
        thisk_kw = M{mat(k_w)}.k;
        thisk_kn = M{mat(k_n)}.k;
        thisk_ks = M{mat(k_s)}.k;
        thisrhoc_p = M{mat(k)}.rhoc_p;
        %
        nonDimFactorx = Deltat/(2.0*thisrhoc_p*Deltax^2);
        nonDimFactory = Deltat/(2.0*thisrhoc_p*Deltay^2);
        % pointer mapping goes row-by-row to assemble Coeff. Matrix
        A(k,k) = 1.0 + (nonDimFactorx*(thisk_ke+2.0*thisk_k+thisk_kw)) + (nonDimFactory*(thisk_kn+2.0*thisk_k+thisk_ks));
        A(k,k_e) = A(k,k_e) - nonDimFactorx*(thisk_ke+thisk_k);
        A(k,k_w) = A(k,k_w) - nonDimFactorx*(thisk_kw+thisk_k);
        A(k,k_n) = A(k,k_n) - nonDimFactory*(thisk_kn+thisk_k);
        A(k,k_s) = A(k,k_s) - nonDimFactory*(thisk_ks+thisk_k);
        %
        Q(k) = (Deltat/(thisrhoc_p*T_r))*q3prime(k);
    end
end

% Center Boundary
i = 1;
j = j_max;
k = pmap(i,j,i_max);
k_e = k + 1;
k_s = k - i_max;
k_n = k_e;
k_w = k_s;
% Treat as normal internal node with special mapping
% matnum = mat(k); % what material
thisk_k = M{mat(k)}.k;
thisk_ke = M{mat(k_e)}.k;
thisk_kw = M{mat(k_w)}.k;
thisk_kn = M{mat(k_n)}.k;
thisk_ks = M{mat(k_s)}.k;
thisrhoc_p = M{mat(k)}.rhoc_p;
%
nonDimFactorx = Deltat/(2.0*thisrhoc_p*Deltax^2);
nonDimFactory = Deltat/(2.0*thisrhoc_p*Deltay^2);
% pointer mapping goes row-by-row to assemble Coeff. Matrix
A(k,k) = 1.0 + (nonDimFactorx*(thisk_ke+2.0*thisk_k+thisk_kw)) + (nonDimFactory*(thisk_kn+2.0*thisk_k+thisk_ks));
A(k,k_e) = A(k,k_e) - nonDimFactorx*(thisk_ke+thisk_k);
A(k,k_s) = A(k,k_s) - nonDimFactory*(thisk_ks+thisk_k);
A(k,k_n) = A(k,k_n) - nonDimFactory*(thisk_kn+thisk_k);
A(k,k_w) = A(k,k_w) - nonDimFactorx*(thisk_kw+thisk_k);
%
Q(k) = (Deltat/(thisrhoc_p*T_r))*q3prime(k);

% Corner Newton Boundary
i = i_max;
j = 1;
k = pmap(i,j,i_max);
k_w = k - 1;
k_n = k + i_max;
% APPLY DOUBLE NEWTON BC 
thisk_k = M{mat(k)}.k;
thish_ke = M{1}.h; % In complete solver, pull in h from the east
thisk_kw = M{mat(k_w)}.k;
thisk_kn = M{mat(k_n)}.k;
thish_ks = M{1}.h; % In complete solver, pull in h from the south
thisrhoc_p = M{mat(k)}.rhoc_p;
%
nonDimFactorx = Deltat/(2.0*thisrhoc_p*Deltax^2);
nonDimFactory = Deltat/(2.0*thisrhoc_p*Deltay^2);
% pointer mapping goes row-by-row to assemble Coeff. Matrix
A(k,k) = 1.0 + (nonDimFactorx*(2.0*thish_ke*Deltax+thisk_k+thisk_kw)) + (nonDimFactory*(thisk_kn+thisk_k+2.0*thish_ks*Deltay));
A(k,k_w) = A(k,k_w) - nonDimFactorx*(thisk_kw+thisk_k);
A(k,k_n) = A(k,k_n) - nonDimFactory*(thisk_kn+thisk_k);
%
Q(k) = (Deltat/(thisrhoc_p*T_r))*q3prime(k) + (thish_ke*Deltat*T_infty/(thisrhoc_p*Deltax*T_r)) + (thish_ks*Deltat*T_infty/(thisrhoc_p*Deltay*T_r));

% Newton Boundaries on Cut Line
%
% Bottom-Left
i = 1;
j = 1;
k = pmap(i,j,i_max);
k_e = k + 1;
k_w = sympmap(i,j,i_max,j_max)-i_max; % kw sym
k_n = k + i_max;
% APPLY NEWTON BC SAME WAY AS IN BOTTOM EDGE
thisk_k = M{mat(k)}.k;
thisk_ke = M{mat(k_e)}.k;
thisk_kw = M{mat(k_w)}.k;
thisk_kn = M{mat(k_n)}.k;
thish_ks = M{1}.h; % In complete solver, pull in h from the south
thisrhoc_p = M{mat(k)}.rhoc_p;
%
nonDimFactorx = Deltat/(2.0*thisrhoc_p*Deltax^2);
nonDimFactory = Deltat/(2.0*thisrhoc_p*Deltay^2);
% pointer mapping goes row-by-row to assemble Coeff. Matrix
A(k,k) = 1.0 + (nonDimFactorx*(thisk_ke+2.0*thisk_k+thisk_kw)) + (nonDimFactory*(thisk_kn+thisk_k+2.0*thish_ks*Deltay));
A(k,k_e) = A(k,k_e) - nonDimFactorx*(thisk_ke+thisk_k);
A(k,k_w) = A(k,k_w) - nonDimFactorx*(thisk_kw+thisk_k);
A(k,k_n) = A(k,k_n) - nonDimFactory*(thisk_kn+thisk_k);
%
Q(k) = (Deltat/(thisrhoc_p*T_r))*q3prime(k) + (thish_ks*Deltat*T_infty/(thisrhoc_p*Deltay*T_r));

% Top-Right
i = i_max;
j = j_max;
k = pmap(i,j,i_max);
k_w = k - 1;
k_n = sympmap(i,j,i_max,j_max)+1; % kn sym
k_s = k - i_max;
% APPLY NEWTON BC SAME WAY AS IN RIGHT EDGE
thisk_k = M{mat(k)}.k;
thish_ke = M{1}.h; % In complete solver, pull in h from the east
thisk_kw = M{mat(k_w)}.k;
thisk_kn = M{mat(k_n)}.k;
thisk_ks = M{mat(k_s)}.k;
thisrhoc_p = M{mat(k)}.rhoc_p;
%
nonDimFactorx = Deltat/(2.0*thisrhoc_p*Deltax^2);
nonDimFactory = Deltat/(2.0*thisrhoc_p*Deltay^2);
% pointer mapping goes row-by-row to assemble Coeff. Matrix
A(k,k) = 1.0 + (nonDimFactorx*(2.0*thish_ke*Deltax+thisk_k+thisk_kw)) + (nonDimFactory*(thisk_kn+2.0*thisk_k+thisk_ks));
A(k,k_w) = A(k,k_w) - nonDimFactorx*(thisk_kw+thisk_k);
A(k,k_n) = A(k,k_n) - nonDimFactory*(thisk_kn+thisk_k);
A(k,k_s) = A(k,k_s) - nonDimFactory*(thisk_ks+thisk_k);
%
Q(k) = (Deltat/(thisrhoc_p*T_r))*q3prime(k) + (thish_ke*Deltat*T_infty/(thisrhoc_p*Deltax*T_r));

%% Solve
% ------------------------------------------------------------------------------

% Init Solution Variables (1D because we use pointer mapping)
T = ones(i_max*j_max,1);
% "Old" Solution for computing residual
T_old = ones(i_max*j_max,1);

% Define variables for time-iteration
residual = 1.0E5; % init residual
epsilon = 1.0E-16; % drive residual down to this value before terminating
myt = 0; % Init time
tTot = 0; % and computation time
p = 0; % track iteration

while (residual > epsilon)
% while (p < 1)
    tStart = tic;

    % Track previous solution vector
    T_old = T;

    % Solve new guess of T
    T = A\(T + Q);

    % Compute the new residual
    residual = norm(T-T_old);
    residual = residual/(i_max*j_max); % Normalize for DOF

    tTot = tTot + toc(tStart);

    % Plot solution
    if mod(p,1000) == 0
        Tplot = T_r*reshape(T, i_max, j_max); % Reshape and Renormalize to Reference Temperature
        figure(1);
        surf(x,y,Tplot);
        ylabel('y');
        xlabel('x');
        title('Temperature Surface');
        drawnow;
    end
    
    myt = myt + Deltat; % Update time (s)
    p = p + 1;

    fprintf(1,'iter = %i, residual = %g\n',p,log10(residual));
end

%% Visualize Results
% ------------------------------------------------------------------------------
Tplot = reshape(T, i_max, j_max);
Tref1 = rot90(Tplot);
Tref1 = Tref1(:,2:end);
Tref1Plot = horzcat(Tplot, Tref1);
Tref2 = rot90(Tref1Plot, 2);
Tref2 = Tref2(1:end-1,:);
Tref2Plot = vertcat(Tref2, Tref1Plot);
Tref2Plot = T_r*Tref2Plot; % Renormalize to Reference Temperature

figure(1);
% Plot flux surface
surf(fullx,fully,Tref2Plot);
ylabel('y [cm]');
xlabel('x [cm]');
zlabel('Temperature [K]')
title('Steady-State Solution for a Uniformly Heated Slab');

%% Store Results
% ------------------------------------------------------------------------------
subsubfolder = [num2str(Deltat),'dt_',num2str(size),'cm'];
plotOut = fullfile(myCWD,subfolder,subsubfolder);
mkdir(plotOut)

% Save Figure
saveas(figure(1),fullfile(plotOut,'tempContour.jpg'));

% And Steady State Solution Matrix
save(fullfile(plotOut,'Tplot.mat'), 'Tref2Plot');

% And Timing Info
tAvg = tTot/p;
fid = fopen(fullfile(plotOut,'time.txt'),'wt');
fprintf(fid, 'Total CPU-time: %s s\nAverage Time per Iteration: %s s\n', string(tTot), string(tAvg));
fclose(fid);

%% Functions
% ------------------------------------------------------------------------------

% Pointer mapping function for 2D->1D transform of solution
function k = pmap(i, j, i_max)
    k = i + (j-1)*i_max;
end

function symk = sympmap(i, j, i_max, j_max)
    symk = (j_max-j+1) + (i_max-i)*j_max;
end