%% 2D Cartesian G-Group Neutron Diffusion
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
UO2 = readtable(fullfile(myCWD,'data\\2Group_PWR.csv'));
H2O = readtable(fullfile(myCWD,'data\\2Group_H2O.csv'));
% UO2 = readtable(fullfile(myCWD,'data\\4Group_10.1GWD_UO2.csv'));
% H2O = readtable(fullfile(myCWD,'data\\4Group_H2O.csv'));

% Cell Array of Materials
M = {UO2, H2O};

% Physical params
totPwr = 3000; % Total Core Power [MW_{th}]
fuelLength = 400; % Total Average Fuel Rod Length [cm]
totLinPwr = 1E6*totPwr/fuelLength; % Total Linear Heat Generation [W/cm]

%% Computational Parameters
% ------------------------------------------------------------------------------
G = 2; % Number of energy groups
% G = 4;

% Define mesh size
% fprintf('Quarter-mesh size') % separate print and input b/c vscode extension
% i_max = input('');
i_max = 65;
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
size = 125; % Domain Size [cm]
maxNodes = i_max*2 - 1; % Total number of nodes in domain
fuelDim = 58; % Fuel Dimensions [Deltax] or [number of nodes]
modDim = ceil((maxNodes-fuelDim)/2);

% Unitless Constants
% Abar = A/L;
% Bbar = B/L;
sizebar = size/L;
% LSig_a = L*Sigma_a;
% LnuSig_f = L*nuSigma_f;

% Calculate step sizes
Deltax = size/(i_max*2 - 1);
Deltay = Deltax;
Deltaxbar = sizebar/(i_max*2 - 1);
Deltaybar = Deltaxbar;

if Deltax > L
    fprintf('Warning: Grid Resolution is Larger than Diffusion Length')
end

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
        if i <= quarterFuelDim && j_max-j+1 <= quarterFuelDim % because of reflection mapping j is flipped
            mat(i,j) = 1; % fuel
            ind(i,j) = pmap(i,j,i_max); % fuel
        else
            mat(i,j) = 2; % moderator
            ind(i,j) = pmap(i,j,i_max); % moderator
        end
    end
end

% Debug Materials Specifications matrix
% matref1 = rot90(mat);
% matref1 = matref1(:,2:end);
% matref1Plot = horzcat(mat, matref1);
% matref2 = rot90(matref1Plot, 2);
% matref2 = matref2(1:end-1,:);
% matref2Plot = vertcat(matref2, matref1Plot);

% File Info
subfolder='results\\project2\\'+string((2*i_max)-1)+'x'+string((2*i_max)-1);
% subfolder='results\\project2\\'+string((2*i_max)-1)+'x'+string((2*i_max)-1)+'_4G';
mkdir(fullfile(myCWD,subfolder));

%% Build Coefficient Matrices
% ------------------------------------------------------------------------------
% Init coeff matrices
H = spalloc(G*i_max*j_max, G*i_max*j_max, 5*G*i_max*j_max); % Sparsely allocate Streaming/Absorption Operator with 5 bands for each energy group
S = spalloc(G*i_max*j_max, G*i_max*j_max, (G-1)*i_max*j_max); % Sparsely allocate Scattering Source Operator with bands for up/downscatter
F = spalloc(G*i_max*j_max, G*i_max*j_max, 1*G*i_max*j_max); % Sparsely allocate Fission Source Operator with 1 band for each energy group

% Init Solution Variables (1D because we use pointer mapping)
phi = ones(G*i_max*j_max,1);
keff = 1;
% "Old" Solution for computing residual
phi_old = ones(G*i_max*j_max,1);
keff_old = 1;

% Define Coefficient Matrix
for group = 1:G
    for i = 2:i_max-1
        for j = 2:j_max-1
            matnum = mat(i, j); % what material
            thisD = M{matnum}.D(group);
            thisSig_R = M{matnum}.Sigma_R(group);

            k = ((group-1)*i_max*j_max) + pmap(i, j, i_max);
            k_e = k + 1;
            k_w = k - 1;
            k_n = k + i_max;
            k_s = k - i_max;

            % pointer mapping goes row-by-row to assemble Coeff. Matrix
            H(k,k) = L*thisSig_R*(1.0) + (thisD/L)*((2.0/Deltaxbar^2) + (2.0/Deltaybar^2));
            H(k,k_e) = (thisD/L)*(-1.0/Deltaxbar^2);
            H(k,k_w) = (thisD/L)*(-1.0/Deltaxbar^2);
            H(k,k_n) = (thisD/L)*(-1.0/Deltaybar^2);
            H(k,k_s) = (thisD/L)*(-1.0/Deltaybar^2);
        end
    end
end
for i = 2:i_max-1
    for j = 2:j_max-1
        matnum = mat(i, j); % what material
        k = pmap(i, j, i_max); % what node
        for group = 1:G % iterate by energy group
            chi_g = M{matnum}.Chi(group); % probability of fission source into this group
            for gprime = 1:G % and by source group 'gprime'

                tok = k + (group-1)*i_max*j_max; 
                fromk = k + (gprime-1)*i_max*j_max;

                nuSig_fgprime = M{matnum}.nuSigma_f(gprime);
                Sig_sgprimeg = M{matnum}.(7+gprime)(group); 

                F(tok,fromk) = L*chi_g*nuSig_fgprime*(1.0);
                S(tok,fromk) = L*Sig_sgprimeg*(1.0);
            end
        end
    end
end
H = sparse(H); % Enforce Coeffs sparse
F = sparse(F);
S = sparse(S);

% Apply BCs
% Left BC
for i = 1:1
    for j = 2:j_max-1 % Avoid corners
        kpos = pmap(i, j, i_max); % what node #
        matnum = mat(i, j); % what material
        for group = 1:G
            thisD = M{matnum}.D(group);
            thisSig_R = M{matnum}.Sigma_R(group);
            chi_g = M{matnum}.Chi(group); % probability of fission source into this group

            k = ((group-1)*i_max*j_max) + kpos;
            k_e = k + 1;
            k_n = k + i_max;
            k_s = k - i_max;
            k_w = ((group-1)*i_max*j_max) + sympmap(i,j,i_max,j_max)-i_max;

            % Treat as normal internal node with special mapping
            H(k,k) = L*thisSig_R*(1.0) + (thisD/L)*((2.0/Deltaxbar^2) + (2.0/Deltaybar^2));
            H(k,k_e) = H(k,k_e) + (thisD/L)*(-1.0/Deltaxbar^2);
            H(k,k_w) = H(k,k_w) + (thisD/L)*(-1.0/Deltaxbar^2);
            H(k,k_n) = H(k,k_n) + (thisD/L)*(-1.0/Deltaybar^2);
            H(k,k_s) = H(k,k_s) + (thisD/L)*(-1.0/Deltaybar^2);
            for gprime = 1:G
                tok = kpos + (group-1)*i_max*j_max; 
                fromk = kpos + (gprime-1)*i_max*j_max;

                nuSig_fgprime = M{matnum}.nuSigma_f(gprime);
                Sig_sgprimeg = M{matnum}.(7+gprime)(group);

                F(tok,fromk) = L*chi_g*nuSig_fgprime*(1.0);
                S(tok,fromk) = L*Sig_sgprimeg*(1.0);                
            end
        end
    end
end
% Right BC
for i = i_max:i_max
    for j = 1:j_max
        kpos = pmap(i, j, i_max); % what node #
        for group = 1:G
            k = ((group-1)*i_max*j_max) + kpos;

            H(k,k) = 1.0;
            % for gprime = 1:G
            %     tok = kpos + (group-1)*i_max*j_max; 
            %     fromk = kpos + (gprime-1)*i_max*j_max;
            %     F(tok,fromk) = 0.0;
            %     S(tok,fromk) = 0.0;                
            % end
        end
    end
end
% Bottom BC
for i = 1:i_max
    for j = 1:1
        kpos = pmap(i, j, i_max); % what node #
        for group = 1:G
            k = ((group-1)*i_max*j_max) + kpos;

            H(k,k) = 1.0;
            % for gprime = 1:G
            %     tok = kpos + (group-1)*i_max*j_max; 
            %     fromk = kpos + (gprime-1)*i_max*j_max;
            %     F(tok,fromk) = 0.0;
            %     S(tok,fromk) = 0.0;                
            % end
        end
    end
end
% Top BC
for i = 2:i_max-1 % Avoid corners
    for j = j_max:j_max
        kpos = pmap(i, j, i_max); % what node #
        matnum = mat(i, j); % what material
        for group = 1:G
            thisD = M{matnum}.D(group);
            thisSig_R = M{matnum}.Sigma_R(group);
            chi_g = M{matnum}.Chi(group); % probability of fission source into this group

            k = ((group-1)*i_max*j_max) + kpos;
            k_e = k + 1;
            k_n = ((group-1)*i_max*j_max) + sympmap(i,j,i_max,j_max)+1;
            k_s = k - i_max;
            k_w = k - 1;

            % Treat as normal internal node with special mapping
            H(k,k) = L*thisSig_R*(1.0) + (thisD/L)*((2.0/Deltaxbar^2) + (2.0/Deltaybar^2));
            H(k,k_e) = H(k,k_e) + (thisD/L)*(-1.0/Deltaxbar^2);
            H(k,k_w) = H(k,k_w) + (thisD/L)*(-1.0/Deltaxbar^2);
            H(k,k_n) = H(k,k_n) + (thisD/L)*(-1.0/Deltaybar^2);
            H(k,k_s) = H(k,k_s) + (thisD/L)*(-1.0/Deltaybar^2);
            for gprime = 1:G
                tok = kpos + (group-1)*i_max*j_max; 
                fromk = kpos + (gprime-1)*i_max*j_max;

                nuSig_fgprime = M{matnum}.nuSigma_f(gprime);
                Sig_sgprimeg = M{matnum}.(7+gprime)(group); 

                F(tok,fromk) = L*chi_g*nuSig_fgprime*(1.0);
                S(tok,fromk) = L*Sig_sgprimeg*(1.0);
            end
        end
    end
end

% Corner/Center Boundary
i = 1;
j = j_max;
kpos = pmap(i,j,i_max);
matnum = mat(i, j); % what material
for group = 1:G
    thisD = M{matnum}.D(group);
    thisSig_R = M{matnum}.Sigma_R(group);
    chi_g = M{matnum}.Chi(group); % probability of fission source into this group

    k = ((group-1)*i_max*j_max) + kpos;
    k_e = k + 1;
    k_s = k - i_max;
    k_n = k_e;
    k_w = k_s;

    % Treat as normal internal node with special mapping
    H(k,k) = L*thisSig_R*(1.0) + (thisD/L)*((2.0/Deltaxbar^2) + (2.0/Deltaybar^2));
    H(k,k_e) = H(k,k_e) + (thisD/L)*(-1.0/Deltaxbar^2);
    H(k,k_w) = H(k,k_w) + (thisD/L)*(-1.0/Deltaxbar^2);
    H(k,k_n) = H(k,k_n) + (thisD/L)*(-1.0/Deltaybar^2);
    H(k,k_s) = H(k,k_s) + (thisD/L)*(-1.0/Deltaybar^2);

    for gprime = 1:G
        tok = kpos + (group-1)*i_max*j_max; 
        fromk = kpos + (gprime-1)*i_max*j_max;

        nuSig_fgprime = M{matnum}.nuSigma_f(gprime);
        Sig_sgprimeg = M{matnum}.(7+gprime)(group); 

        F(tok,fromk) = L*chi_g*nuSig_fgprime*(1.0);
        S(tok,fromk) = L*Sig_sgprimeg*(1.0);
    end
end

%% Power-Iteration Script
% ------------------------------------------------------------------------------
% Compute evolution operator initially to minimize work in loop
% Amat = inv(H-S)*F % Slowww
% Amat = (H-S)\F; % Memory-Intensive
Amat = (H-S)\F; % Enforce Memory Clearing (since MATLAB is weird about it)
fprintf('bruh')
pause()
clearvars H S F

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

    % Solve new guess of k (optional)
    % phiT = transpose(phi);
    % keff = phiT * (Amat * phi); % Search Dominant Eigenvalue

    % Compute the new residual
    residual = norm(phi-phi_old);
    residual = residual/(i_max*j_max); % Normalize for DOF

    % Plot solution
    if mod(iter,10) == 0
        mygroup = G; % Look at slowest group because characteristic features
        slowPlot = reshape(phi(1+(mygroup-1)*i_max*j_max:mygroup*i_max*j_max), i_max, j_max);
        figure(1);
        % Plot flux surface
        surf(x,y,slowPlot);
        ylabel('y');
        xlabel('x');
        title('Thermal Flux Surface');
        drawnow;
    end

    tTot = tTot + toc(tStart);

    fprintf(1,'iter = %i, residual = %g\n',iter,log10(residual));
    iter = iter + 1;
end

% Final Value of k
phiT = transpose(phi);
keff = phiT * (Amat * phi); % Search Dominant Eigenvalue

%% Visualize Results
% ------------------------------------------------------------------------------
phiOrganized = zeros(G, i_max*j_max);
for mygroup=1:G
    phiOrganized(mygroup,:) = phi(1+(mygroup-1)*i_max*j_max:mygroup*i_max*j_max);
end

phiPlotOrganized = zeros(G, i_max*2 - 1, j_max*2 - 1);
for mygroup = 1:G
    phiPlot = reshape(phiOrganized(mygroup,:), i_max, j_max);
    phiref1 = rot90(phiPlot);
    phiref1 = phiref1(:,2:end);
    phiref1Plot = horzcat(phiPlot, phiref1);
    phiref2 = rot90(phiref1Plot, 2);
    phiref2 = phiref2(1:end-1,:);
    phiref2Plot = vertcat(phiref2, phiref1Plot);

    phiPlotOrganized(mygroup, :, :) = phiref2Plot;
end
phiPlotOrganized = phiPlotOrganized/(sum(phiPlotOrganized,'all')*Deltax*Deltay); % Renormalize by true area under the surface

% Physical Units
phiPlotOrganized = phiPlotOrganized*totLinPwr; % Normalize to units of [W/cc]

for mygroup = 1:G
    figure(mygroup);
    % Plot flux surface
    phiPlot = reshape(phiPlotOrganized(mygroup, :, :), i_max*2 - 1, j_max*2 - 1);
    surf(fullx,fully,phiPlot);
    ylabel('y [cm]');
    xlabel('x [cm]');
    zlabel('Volumetric Heat Rate [W/cc]')
    title(['Group ',num2str(mygroup),' Flux Surface']);
end

fprintf(1,'keff = %f\n',keff);

%% Store Results
% ------------------------------------------------------------------------------

% Redimensionalize for data output
fuelSize = fuelDim*Deltax; 
modThick = modDim*Deltax;

% Ensure we have outpit directories
plotsfolder = [num2str(fuelDim),'_',num2str(modDim),'_',num2str(size)];
plotOut = fullfile(myCWD,subfolder,plotsfolder);
mkdir(plotOut)
resultmat = 'myTable.mat';
resultOut = fullfile(myCWD,subfolder,resultmat);

newData = table(fuelSize, modThick, size, keff,'VariableNames', {'fuelSize', 'modThick', 'size', 'keff'});

if exist(resultOut, 'file') == 2
    % Load the existing table
    existingTable = load(resultOut);
    
    % Assuming the table is stored in a variable within the .mat file with the same name as the file
    % Adjust this if your table variable name is different
    existingTableVarName = genvarname(strrep(resultOut, '.mat', ''));
    
    % Append the new data
    myTable = [existingTable.myTable; newData];
    
    % Save the updated table back to the file
    save(resultOut, 'myTable'); % Save the updated table under a new name or overwrite
else
    % If the table file does not exist, create a new one
    myTable = newData;
    save(resultOut, 'myTable');
end

fprintf(1, 'Slab Dimension: %f cm\n', fuelSize);
fprintf(1, 'Moderator Thickness: %f cm\n', modThick);

% Save Plots for Each Energy Group
for mygroup = 1:G
    saveas(figure(mygroup),fullfile(plotOut,['G',num2str(mygroup),'FluxContour.jpg']));
end

% And Solution Matrix
save(fullfile(plotOut,'phiPlot.mat'), 'phiPlotOrganized');

% And Timing Info
tAvg = tTot/iter;
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