%% 2D Coupled SMR
% ------------------------------------------------------------------------------
% Adam Buchalter
%
% desc
% ------------------------------------------------------------------------------
clear all; close all; clc;
myCWD = pwd;

%% User Input
% ------------------------------------------------------------------------------

% Thermal Data
FUL = readtable(fullfile(myCWD,'proj6Data\\thermalData\\Fuel_Lump.csv'));
FU2 = readtable(fullfile(myCWD,'proj6Data\\thermalData\\Fuel_Lump.csv'));
MOD = readtable(fullfile(myCWD,'proj6Data\\thermalData\\H2O_Lump.csv'));
% Cell Array of Materials
%     1    2    3    4
M = {FUL, FU2, MOD};
MODindeces = [3]; % Which Indeces in the Cell Array are Coolant?

% Neutronics Interpolation Tables in Same Order as Materials Array
myLibraries = {
    fullfile(myCWD,'proj6Data\\neutronData\\interpTables.mat'), % FUL 
    fullfile(myCWD,'proj6Data\\neutronData\\interpTables2.mat'), % FU2 
    fullfile(myCWD,'proj6Data\\neutronData\\interpTablesMod.mat') % MOD
};

% Quarter-Core Layout
LAYOUT = [
    2 2 2 1 3
    2 2 1 1 3
    2 1 1 3 3
    1 1 3 3 3
    3 3 3 3 3
];

TF_nom = 850.0; % Nominal Fuel Temperature [K]
TM_nom = 587.0; % Nominal Moderator Temperature [K]

myBORHi = 2100; % Boron Concentration [ppm] that yields k-eff > 1
myBORLo = 2230; % k-eff < 1
% myBOR = 4780; % Soluble Poison Concentration

% File Info
subfolder='results\\project6\\shellFuel_TM587';
mkdir(fullfile(myCWD,subfolder));

%% Import Data
% ------------------------------------------------------------------------------

% Neutron Data
shapeMats = size(M); % Number of Materials

% Init Fitting Tables as Cell Arrays to Generalize
global TFs;
TFs=cell(shapeMats);
global TMs;
TMs=cell(shapeMats);
global BORs;
BORs=cell(shapeMats);
global D1;
D1=cell(shapeMats);
global D2;
D2=cell(shapeMats);
global Sigma_R1;
Sigma_R1=cell(shapeMats);
global Sigma_R2;
Sigma_R2=cell(shapeMats);
global Sigma_a1;
Sigma_a1=cell(shapeMats);
global Sigma_a2;
Sigma_a2=cell(shapeMats);
% global Sigma_s11;
% Sigma_s11=cell(shapeMats);
global Sigma_s12;
Sigma_s12=cell(shapeMats);
global Sigma_s21;
Sigma_s21=cell(shapeMats);
% global Sigma_s22;
% Sigma_s22=cell(shapeMats);
global nuSigma_f1;
nuSigma_f1=cell(shapeMats);
global nuSigma_f2;
nuSigma_f2=cell(shapeMats);

for index=1:shapeMats(2)
    myMat = load(myLibraries{index});
    TFs{index} = myMat.TFs;
    TMs{index} = myMat.TMs;
    BORs{index} = myMat.BORs;
    D1{index} = myMat.D1;
    D2{index} = myMat.D2;
    Sigma_R1{index} = myMat.Sigma_R1;
    Sigma_R2{index} = myMat.Sigma_R2;
    Sigma_a1{index} = myMat.Sigma_a1;
    Sigma_a2{index} = myMat.Sigma_a2;
    % Sigma_s11{index} = myMat.Sigma_s11;
    Sigma_s12{index} = myMat.Sigma_s12;
    Sigma_s21{index} = myMat.Sigma_s21;
    % Sigma_s22{index} = myMat.Sigma_s22;
    nuSigma_f1{index} = myMat.nuSigma_f1;
    nuSigma_f2{index} = myMat.nuSigma_f2;
end
CHI = [1 0]; % Assume ALL Neutrons Birthed in Fast Group

% Physical params
totPwr = 250*1E6; % Thermal Output [MW_th -> W_th]
% NbyN = 17; % 17x17 Fuel Assembly
fuelLength = 78.74*2.54; % Active Height [in -> cm]
pitch_Assy = 8.466*2.54; % Assembly Pitch [in -> cm]
totLinPwr = totPwr/fuelLength; % Total Linear Heat Generation [W/cm]
% myVel = 2.7*12*2.54; % Average in-core flow velocity [ft/s -> cm/s]
myVel = 2.1*12*2.54;
T_in = (497-32)*(5/9) + 273.15; % Inlet Temperature [F->K]
T_out = (597-32)*(5/9) + 273.15; % Outlet Temperature [F->K]

% Area Corrections
N_FP = 264; % Fuel Pins per Assy
N_GT = 25; % Guide and Instrument Tubes per Assy
OD_GT = 0.482*2.54; % Guide Tube Outer Diameter [in -> cm]
pitch_Assy = 21.50; % Assembly Pitch [cm]
c = 0.024*2.54; % Clad Thickness [in -> cm]
g = 0.0065*2.54; % Pellet-Clad Gap [in -> cm]
s = 0.3195*2.54; % Fuel Pellet Diameter [in -> cm]
R = s/2; % Fuel-pellet Radius [cm]
assyArea = pitch_Assy^2; % Fuel Assembly Area [cm^2]
fuelArea = N_FP*pi*R^2; % Area Occupied by Fuel [cm^2]
pinsArea = N_FP*pi*(R+g+c)^2 + N_GT*pi*(OD_GT/2)^2; % Area Occupied by not Water [cm^2]
waterArea = assyArea-pinsArea; % Area Occupied by Water [cm^2]
fuelCorr = fuelArea/assyArea; % Fraction of Assembly Occupied by Fuel
modCorr = waterArea/assyArea; % Fraction of Assembly Occupied by Water

TFUL = TF_nom;
TMOD = TM_nom;

%% Intelligently Find Smallest Diffusion Length for Resolution-Setting
% ------------------------------------------------------------------------------
% The smallest diffusion length is important because if our cells are larger
% than the diffusion length for a particular group then we cannot reliably
% capture transport phenomena

L = 50; % [cm] Arbitrarily Large Initial Guess
for m = 1:shapeMats(2)
    for i = 1:size(TFs{m})
        for j = 1:size(TMs{m})
            for k = 1:size(BORs{m})
                thisL = sqrt(D1{m}(i,j,k)/Sigma_a1{m}(i,j,k)); % Characteristic Diffusion Length [cm]
                if thisL < L % If smaller, use as new normalization length
                    L = thisL;
                end
                thisL = sqrt(D2{m}(i,j,k)/Sigma_a2{m}(i,j,k)); % Characteristic Diffusion Length [cm]
                if thisL < L % If smaller, use as new normalization length
                    L = thisL;
                end
            end
        end
    end
end

fprintf('Minimum Neutron Diffusion Length: %.2f cm\n',L);

Nmin = ceil(pitch_Assy/L);
if mod(Nmin, 2) == 0
    Nmin = Nmin+1;
end

Deltax = pitch_Assy/Nmin;
Deltay = Deltax;

fprintf('Spatial Resolution: %.2f cm\n',Deltax);

fprintf('\nAuto-Scaling Domain...\n');
sizeLayout = size(LAYOUT);
iArrMax = sizeLayout(1);
jArrMax = sizeLayout(2);
i_max = (iArrMax-1)*Nmin + (Nmin+1)/2;
j_max = (jArrMax-1)*Nmin + (Nmin+1)/2;

Domain = zeros(i_max, j_max);
halfSize = (Nmin+1)/2;
% Inside Corner
iArr = 1;
jArr = 1;
Domain(1:halfSize,1:halfSize) = LAYOUT(iArr, jArr);
% Branch Cut 1
jArr = 1;
for iArr = 2:iArrMax
    Domain(halfSize+(iArr-2)*Nmin+1:halfSize+(iArr-1)*Nmin, 1:halfSize) = LAYOUT(iArr, jArr);
end
% Branch Cut 2
iArr = 1;
for jArr = 2:jArrMax
    Domain(1:halfSize, halfSize+(jArr-2)*Nmin+1:halfSize+(jArr-1)*Nmin) = LAYOUT(iArr, jArr);
end
% Full Assemblies
for iArr = 2:iArrMax
    for jArr = 2:jArrMax
        Domain(halfSize+(iArr-2)*Nmin+1:halfSize+(iArr-1)*Nmin, halfSize+(jArr-2)*Nmin+1:halfSize+(jArr-1)*Nmin) = LAYOUT(iArr, jArr);
    end
end

for i = 1:i_max
    for j = 1:j_max
        k = pmap(i, j, j_max);
        mat(k) = Domain(i, j);
    end
end

%% Computational Parameters
% ------------------------------------------------------------------------------
% Bulk Convective Fluid Temperature
T_infty = TM_nom; % [K]

G = 2; % Number of energy groups

%% Nondimensional Domain Prep
% ------------------------------------------------------------------------------

% Define physical domain
size = Deltax*(i_max-1)*2; % FULL Domain Size [cm]

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

%% First, Solve for Fluid Flow (Invariate)
% ------------------------------------------------------------------------------
fprintf('\n---Solving for Coolant Flow Profile---\n');
fprintf('Initializing...\n');
% Init coeff matrices
A = spalloc(i_max*j_max, i_max*j_max, 5*i_max*j_max); % Sparsely allocate Diffusion Operator with 5 bands
S = zeros(i_max*j_max, 1); % Allocate Source Term for Dirichlet Conditions

% Pre-calc these
factorx = 1/Deltax^2;
factory = 1/Deltay^2;

% Define Matrices
fprintf('Building Matrices...\n');
for i = 2:i_max-1
    for j = 2:j_max-1
        k = pmap(i, j, j_max);
        k_e = k + 1;
        k_w = k - 1;
        k_n = k - j_max;
        k_s = k + j_max;

        thisDomain = Domain(i, j);
        if any(MODindeces == mat(k)) % If coolant
            A(k,k) = -2.0*(factorx + factory);
            A(k,k_e) = factorx;
            A(k,k_w) = factorx;
            A(k,k_n) = factory;
            A(k,k_s) = factory;
        else % Fuel
            A(k, k) = 1.0;
            S(k) = myVel; % Bulk Core Velocity Source [cm/s]
        end
    end
end

% Apply BCs
% Left BC
for i = 2:i_max-1 % Avoid corners
    for j = 1:1
        k = pmap(i,j,j_max);
        k_e = k + 1;
        k_n = k - j_max;
        k_s = k + j_max;
        k_w = sympmap(i,j,j_max)+j_max;

        thisDomain = Domain(i, j);
        if any(MODindeces == mat(k)) % If coolant
            A(k,k) = -2.0*(factorx + factory);
            A(k,k_e) = A(k,k_e) + factorx;
            A(k,k_w) = A(k,k_w) + factorx;
            A(k,k_n) = A(k,k_n) + factory;
            A(k,k_s) = A(k,k_s) + factory;
        else % Fuel
            A(k, k) = 1.0;
            S(k) = myVel; % Bulk Core Velocity Source [cm/s]
        end
    end
end
% Right BC -- No Slip
for i = 1:i_max
    for j = j_max:j_max
        k = pmap(i,j,j_max);
        A(k, k) = 1.0;
        S(k) = 0.0;
    end
end
% Bottom BC -- No Slip
for i = i_max:i_max
    for j = 1:j_max
        k = pmap(i,j,j_max);
        A(k, k) = 1.0;
        S(k) = 0.0;
    end
end
% Top BC
for i = 1:1
    for j = 2:j_max-1 % Avoid corners
        k = pmap(i,j,j_max);
        k_e = k + 1;
        k_w = k - 1;
        k_s = k + j_max;
        k_n = sympmap(i,j,j_max)+1;

        thisDomain = Domain(i, j);
        if any(MODindeces == mat(k)) % If coolant
            A(k,k) = -2.0*(factorx + factory);
            A(k,k_e) = A(k,k_e) + factorx;
            A(k,k_w) = A(k,k_w) + factorx;
            A(k,k_n) = A(k,k_n) + factory;
            A(k,k_s) = A(k,k_s) + factory;
        else % Fuel
            A(k, k) = 1.0;
            S(k) = myVel; % Bulk Core Velocity Source [cm/s]
        end
    end
end

% Center Boundary
i = 1;
j = 1;
k = pmap(i,j,j_max);
k_e = k + 1;
k_s = k + i_max;
k_n = k_e;
k_w = k_s;

thisDomain = Domain(i, j);
if any(MODindeces == mat(k)) % If coolant
    A(k,k) = -2.0*(factorx + factory);
    A(k,k_e) = A(k,k_e) + factorx;
    A(k,k_w) = A(k,k_w) + factorx;
    A(k,k_n) = A(k,k_n) + factory;
    A(k,k_s) = A(k,k_s) + factory;
else % Fuel
    A(k, k) = 1.0;
    S(k) = myVel; % Bulk Core Velocity Source [cm/s]
end

fprintf('Solving for %i degrees of freedom...\n', i_max*j_max);
w = A\S; % Axial Flow Field [cm/s] -> forms a heat removal term
fprintf('Complete!\n');
% wPlot = reshape(w, i_max, j_max);
% figure(1);
% surf(x, y, wPlot)
% ylabel('y');
% xlabel('x');
% title('Flow Surface');
% drawnow;

% pause()

clear A S

%% Iterative Solver
% ------------------------------------------------------------------------------
fprintf('\n===ITERATIVE SOLVER===\n');
fprintf('Initializing...\n');

% Init Solution Variables (1D because we use pointer mapping)
T = zeros(i_max*j_max,1);
q3prime = zeros(i_max*j_max,1);
phi = ones(G*i_max*j_max,1);
keff = 1;
% "Old" Solution for computing residual
phi_old = ones(G*i_max*j_max,1);
keff_old = 1;

% Nondim
Deltaxbar = Deltax/L;
Deltaybar = Deltay/L;

% Inital Value Sets
for i = 1:i_max
    for j = 1:j_max
        thisDomain = Domain(i, j);
        k = pmap(i, j, j_max); % what node #
        if thisDomain == 1 % Fuel
            T(k) = TF_nom;
        else % if thisDomain == 2 % Moderator
            T(k) = TM_nom;
        end
    end
end

keff_iter = 2.0;
kEpsilon = 0.00001;

while abs(keff_iter - 1.0) > kEpsilon
% for bigBigIter = 1:10
    myBOR = (myBORHi + myBORLo)/2.0;

    iter=2;

    while (iter>1)
    % for bigIter = 1:10
        %% Iterative Solver: Neutronics
        % ------------------------------------------------------------------------------

        fprintf('\n---Solving for Neutron Flux Profile---\n');

        % Define Matrices
        fprintf('Building Matrices...\n');

        % Init Neutronics Matrices
        H = spalloc(G*i_max*j_max, G*i_max*j_max, 5*G*i_max*j_max); % Sparsely allocate Streaming/Absorption Operator with 5 bands for each energy group
        S = spalloc(G*i_max*j_max, G*i_max*j_max, (G-1)*i_max*j_max); % Sparsely allocate Scattering Source Operator with bands for up/downscatter
        F = spalloc(G*i_max*j_max, G*i_max*j_max, 1*G*i_max*j_max); % Sparsely allocate Fission Source Operator with 1 band for each energy group

        for i = 2:i_max-1
            for j = 2:j_max-1
                kpos = pmap(i, j, j_max); % what node #
                for group = 1:G
                    k = ((group-1)*i_max*j_max) + kpos;
                    k_e = k + 1;
                    k_w = k - 1;
                    k_n = k - j_max;
                    k_s = k + j_max;

                    thisMat = Domain(i, j);
                    thisT = T(kpos);

                    thisD = myD(thisT, TMOD, myBOR, group, thisMat);
                    thisSig_R = mySigma_R(thisT, TMOD, myBOR, group, thisMat);

                    % pointer mapping goes row-by-row to assemble Coeff. Matrix
                    H(k,k) = L*thisSig_R*(1.0) + (thisD/L)*((2.0/Deltaxbar^2) + (2.0/Deltaybar^2));
                    H(k,k_e) = H(k,k_e) + (thisD/L)*(-1.0/Deltaxbar^2);
                    H(k,k_w) = H(k,k_w) + (thisD/L)*(-1.0/Deltaxbar^2);
                    H(k,k_n) = H(k,k_n) + (thisD/L)*(-1.0/Deltaybar^2);
                    H(k,k_s) = H(k,k_s) + (thisD/L)*(-1.0/Deltaybar^2);

                    chi_g = CHI(group);
                    for gprime = 1:G
                        % tok = kpos + (group-1)*i_max*j_max;
                        tok = k;
                        fromk = kpos + (gprime-1)*i_max*j_max;

                        nuSig_fgprime = mynuSigma_f(thisT, TMOD, myBOR, gprime, thisMat);
                        Sig_sgprimeg = mySigma_s(thisT, TMOD, myBOR, gprime, group, thisMat);

                        F(tok,fromk) = L*chi_g*nuSig_fgprime*(1.0);
                        S(tok,fromk) = L*Sig_sgprimeg*(1.0);

                    end
                end
            end
        end

        % Apply BCs
        % Left BC
        for i = 2:i_max-1 % Avoid corners
            for j = 1:1
                kpos = pmap(i, j, j_max); % what node #
                kpossym = sympmap(i,j,j_max);
                for group = 1:G
                    k = ((group-1)*i_max*j_max) + kpos;
                    k_e = k + 1;
                    k_n = k - j_max;
                    k_s = k + j_max;
                    k_w = ((group-1)*i_max*j_max) + kpossym + j_max;

                    thisMat = Domain(i, j);
                    thisT = T(kpos);

                    thisD = myD(thisT, TMOD, myBOR, group, thisMat);
                    thisSig_R = mySigma_R(thisT, TMOD, myBOR, group, thisMat);

                    % pointer mapping goes row-by-row to assemble Coeff. Matrix
                    H(k,k) = L*thisSig_R*(1.0) + (thisD/L)*((2.0/Deltaxbar^2) + (2.0/Deltaybar^2));
                    H(k,k_e) = H(k,k_e) + (thisD/L)*(-1.0/Deltaxbar^2);
                    H(k,k_w) = H(k,k_w) + (thisD/L)*(-1.0/Deltaxbar^2);
                    H(k,k_n) = H(k,k_n) + (thisD/L)*(-1.0/Deltaybar^2);
                    H(k,k_s) = H(k,k_s) + (thisD/L)*(-1.0/Deltaybar^2);

                    chi_g = CHI(group);
                    for gprime = 1:G
                        % tok = kpos + (group-1)*i_max*j_max;
                        tok = k;
                        fromk = kpos + (gprime-1)*i_max*j_max;

                        nuSig_fgprime = mynuSigma_f(thisT, TMOD, myBOR, gprime, thisMat);
                        Sig_sgprimeg = mySigma_s(thisT, TMOD, myBOR, gprime, group, thisMat);

                        F(tok,fromk) = L*chi_g*nuSig_fgprime*(1.0);
                        S(tok,fromk) = L*Sig_sgprimeg*(1.0);

                    end
                end
            end
        end
        % Right BC
        for i = 1:i_max
            for j = j_max:j_max
                kpos = pmap(i, j, j_max); % what node #
                for group = 1:G
                    k = ((group-1)*i_max*j_max) + kpos;
                    
                    H(k,k) = 1.0;
                end
            end
        end
        % Bottom BC
        for i = i_max:i_max
            for j = 1:j_max
                kpos = pmap(i, j, j_max); % what node #
                for group = 1:G
                    k = ((group-1)*i_max*j_max) + kpos;

                    H(k,k) = 1.0;
                end
            end
        end
        % Top BC
        for i = 1:1
            for j = 2:j_max-1 % Avoid corners
                kpos = pmap(i, j, j_max); % what node #
                kpossym = sympmap(i,j,j_max);
                for group = 1:G
                    k = ((group-1)*i_max*j_max) + kpos;
                    k_e = k + 1;
                    k_w = k - 1;
                    k_s = k + j_max;
                    k_n = ((group-1)*i_max*j_max) + kpossym + 1;

                    thisMat = Domain(i, j);
                    thisT = T(kpos);

                    thisD = myD(thisT, TMOD, myBOR, group, thisMat);
                    thisSig_R = mySigma_R(thisT, TMOD, myBOR, group, thisMat);

                    % pointer mapping goes row-by-row to assemble Coeff. Matrix
                    H(k,k) = L*thisSig_R*(1.0) + (thisD/L)*((2.0/Deltaxbar^2) + (2.0/Deltaybar^2));
                    H(k,k_e) = H(k,k_e) + (thisD/L)*(-1.0/Deltaxbar^2);
                    H(k,k_w) = H(k,k_w) + (thisD/L)*(-1.0/Deltaxbar^2);
                    H(k,k_n) = H(k,k_n) + (thisD/L)*(-1.0/Deltaybar^2);
                    H(k,k_s) = H(k,k_s) + (thisD/L)*(-1.0/Deltaybar^2);

                    chi_g = CHI(group);
                    for gprime = 1:G
                        % tok = kpos + (group-1)*i_max*j_max;
                        tok = k;
                        fromk = kpos + (gprime-1)*i_max*j_max;

                        nuSig_fgprime = mynuSigma_f(thisT, TMOD, myBOR, gprime, thisMat);
                        Sig_sgprimeg = mySigma_s(thisT, TMOD, myBOR, gprime, group, thisMat);

                        F(tok,fromk) = L*chi_g*nuSig_fgprime*(1.0);
                        S(tok,fromk) = L*Sig_sgprimeg*(1.0);

                    end
                end
            end
        end

        % Center Boundary
        i = 1;
        j = 1;
        kpos = pmap(i, j, j_max); % what node #
        for group = 1:G
            k = ((group-1)*i_max*j_max) + kpos;
            k_e = k + 1;
            k_s = k + i_max;
            k_n = k_e;
            k_w = k_s;

            thisMat = Domain(i, j);
            thisT = T(kpos);

            thisD = myD(thisT, TMOD, myBOR, group, thisMat);
            thisSig_R = mySigma_R(thisT, TMOD, myBOR, group, thisMat);

            % pointer mapping goes row-by-row to assemble Coeff. Matrix
            H(k,k) = L*thisSig_R*(1.0) + (thisD/L)*((2.0/Deltaxbar^2) + (2.0/Deltaybar^2));
            H(k,k_e) = H(k,k_e) + (thisD/L)*(-1.0/Deltaxbar^2);
            H(k,k_w) = H(k,k_w) + (thisD/L)*(-1.0/Deltaxbar^2);
            H(k,k_n) = H(k,k_n) + (thisD/L)*(-1.0/Deltaybar^2);
            H(k,k_s) = H(k,k_s) + (thisD/L)*(-1.0/Deltaybar^2);

            chi_g = CHI(group);
            for gprime = 1:G
                % tok = kpos + (group-1)*i_max*j_max;
                tok = k;
                fromk = kpos + (gprime-1)*i_max*j_max;

                nuSig_fgprime = mynuSigma_f(thisT, TMOD, myBOR, gprime, thisMat);
                Sig_sgprimeg = mySigma_s(thisT, TMOD, myBOR, gprime, group, thisMat);

                F(tok,fromk) = L*chi_g*nuSig_fgprime*(1.0);
                S(tok,fromk) = L*Sig_sgprimeg*(1.0);

            end
        end

        %% Iterative Solver: Neutronics Power-Iteration
        % ------------------------------------------------------------------------------
        % Compute evolution operator initially to minimize work in loop
        % Amat = inv(H-S)*F % Slowww
        % Amat = (H-S)\F; % Memory-Intensive
        fprintf('Unifying Matrices...\n');
        Amat = (H-S)\F; % Enforce Memory Clearing (since MATLAB is weird about it)
        clearvars H S F

        % Init iteration vars
        tTot = 0;
        iter = 0;
        residual = 1.0E5; % init residual
        epsilon = 1.0E-16; % drive residual down to this value before terminating

        % Begin iteration
        fprintf('Solving for %i degrees of freedom...\n', G*i_max*j_max);
        while (residual > epsilon)
        % while (iter<100)

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
            if mod(iter,100) == 0
                mygroup = G; % Look at slowest group because characteristic features
                slowPlot = reshape(phi(1+(mygroup-1)*i_max*j_max:mygroup*i_max*j_max), i_max, j_max);
                % slowPlot = reshape(phi(1:i_max*j_max), i_max, j_max) + reshape(phi(1+i_max*j_max:2*i_max*j_max), i_max, j_max);
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
        fprintf('Complete!\n');

        % Final Value of k
        phiT = transpose(phi);
        keff = phiT * (Amat * phi); % Search Dominant Eigenvalue
        fprintf('k_eff = %.5f \n',keff);
        fprintf('Boron Conc = %.3f ppm\n',myBOR);

        q3prime = phi(1:i_max*j_max)+phi(1+i_max*j_max:2*i_max*j_max);
        q3prime = (totLinPwr/4)*q3prime/(sum(q3prime,'all')*Deltax*Deltay);

        clearvars Amat

        % pause()

        %% Iterative Solver: Heat Diffusion
        % ------------------------------------------------------------------------------

        fprintf('\n---Solving for Temperature Profile---\n');
        MOD_rhoc_p = M{2}.rhoCp;
        dTdz = (T_out-T_in)/fuelLength;
        % h = 0.05;
        h = 3.0;

        % Define Matrices
        fprintf('Building Matrices...\n');
        % Init Thermal Matrices
        A = spalloc(i_max*j_max, i_max*j_max, 5*i_max*j_max); % Sparsely allocate Diffusion Operator with 5 bands
        Q = zeros(i_max*j_max, 1); % Allocate Heat Source/Removal Vector
        
        for i = 2:i_max-1
            for j = 2:j_max-1
                k = pmap(i, j, j_max); % what node #
                k_e = k + 1;
                k_w = k - 1;
                k_n = k - j_max;
                k_s = k + j_max;

                thisk_k = M{mat(k)}.k;
                thisk_ke = M{mat(k_e)}.k;
                thisk_kw = M{mat(k_w)}.k;
                thisk_kn = M{mat(k_n)}.k;
                thisk_ks = M{mat(k_s)}.k;
                % thisrhoc_p = M{mat(k)}.rhoCp;

                if any(MODindeces == mat(k)) % If coolant
                    heatremoval = MOD_rhoc_p*dTdz*w(k);
                    if mat(k_e) == 1
                        flux_e = h/Deltax;
                    else
                        flux_e = (thisk_ke+thisk_k)/Deltax^2;
                    end
                    if mat(k_w) == 1
                        flux_w = h/Deltax;
                    else
                        flux_w = (thisk_kw+thisk_k)/Deltax^2;
                    end
                    if mat(k_n) == 1
                        flux_n = h/Deltay;
                    else
                        flux_n = (thisk_kn+thisk_k)/Deltay^2;
                    end
                    if mat(k_s) == 1
                        flux_s = h/Deltay;
                    else
                        flux_s = (thisk_ks+thisk_k)/Deltay^2;
                    end
                else % Fuel
                    heatremoval = modCorr*MOD_rhoc_p*dTdz*w(k);
                    if any(MODindeces == mat(k_e))
                        flux_e = h/Deltax;
                    else
                        flux_e = (thisk_ke+thisk_k)/Deltax^2;
                    end
                    if any(MODindeces == mat(k_w))
                        flux_w = h/Deltax;
                    else
                        flux_w = (thisk_kw+thisk_k)/Deltax^2;
                    end
                    if any(MODindeces == mat(k_n))
                        flux_n = h/Deltay;
                    else
                        flux_n = (thisk_kn+thisk_k)/Deltay^2;
                    end
                    if any(MODindeces == mat(k_s))
                        flux_s = h/Deltay;
                    else
                        flux_s = (thisk_ks+thisk_k)/Deltay^2;
                    end
                end
                A(k,k) = flux_e + flux_w + flux_n + flux_s;
                A(k,k_e) = A(k,k_e) - flux_e;
                A(k,k_w) = A(k,k_w) - flux_w;
                A(k,k_n) = A(k,k_n) - flux_n;
                A(k,k_s) = A(k,k_s) - flux_s;
                Q(k) = (q3prime(k) - heatremoval);
            end
        end

        % Apply BCs
        % Left BC
        for i = 2:i_max-1 % Avoid corners
            for j = 1:1
                k = pmap(i, j, j_max); % what node #
                ksym = sympmap(i,j,j_max);
                k_e = k + 1;
                k_n = k - j_max;
                k_s = k + j_max;
                k_w = ksym + j_max;

                % thisMat = Domain(k); % what material
                % matnum = mat(k); % what material
                thisk_k = M{mat(k)}.k;
                thisk_ke = M{mat(k_e)}.k;
                thisk_kw = M{mat(k_w)}.k;
                thisk_kn = M{mat(k_n)}.k;
                thisk_ks = M{mat(k_s)}.k;
                % thisrhoc_p = M{mat(k)}.rhoCp;

                if any(MODindeces == mat(k)) % If coolant
                    heatremoval = MOD_rhoc_p*dTdz*w(k);
                    if mat(k_e) == 1
                        flux_e = h/Deltax;
                    else
                        flux_e = (thisk_ke+thisk_k)/Deltax^2;
                    end
                    if mat(k_w) == 1
                        flux_w = h/Deltax;
                    else
                        flux_w = (thisk_kw+thisk_k)/Deltax^2;
                    end
                    if mat(k_n) == 1
                        flux_n = h/Deltay;
                    else
                        flux_n = (thisk_kn+thisk_k)/Deltay^2;
                    end
                    if mat(k_s) == 1
                        flux_s = h/Deltay;
                    else
                        flux_s = (thisk_ks+thisk_k)/Deltay^2;
                    end
                else % Fuel
                    heatremoval = modCorr*MOD_rhoc_p*dTdz*w(k);
                    if any(MODindeces == mat(k_e))
                        flux_e = h/Deltax;
                    else
                        flux_e = (thisk_ke+thisk_k)/Deltax^2;
                    end
                    if any(MODindeces == mat(k_w))
                        flux_w = h/Deltax;
                    else
                        flux_w = (thisk_kw+thisk_k)/Deltax^2;
                    end
                    if any(MODindeces == mat(k_n))
                        flux_n = h/Deltay;
                    else
                        flux_n = (thisk_kn+thisk_k)/Deltay^2;
                    end
                    if any(MODindeces == mat(k_s))
                        flux_s = h/Deltay;
                    else
                        flux_s = (thisk_ks+thisk_k)/Deltay^2;
                    end
                end
                A(k,k) = flux_e + flux_w + flux_n + flux_s;
                A(k,k_e) = A(k,k_e) - flux_e;
                A(k,k_w) = A(k,k_w) - flux_w;
                A(k,k_n) = A(k,k_n) - flux_n;
                A(k,k_s) = A(k,k_s) - flux_s;
                Q(k) = (q3prime(k) - heatremoval);
            end
        end
        % Right BC
        for i = 1:i_max
            for j = j_max:j_max
                k = pmap(i, j, j_max); % what node #
                A(k,k) = 1.0;
                Q(k) = T_infty;
            end
        end
        % Bottom BC
        for i = i_max:i_max
            for j = 1:j_max
                k = pmap(i, j, j_max); % what node #
                A(k,k) = 1.0;
                Q(k) = T_infty;
            end
        end

        % Top BC
        for i = 1:1
            for j = 2:j_max-1 % Avoid corners
                k = pmap(i, j, j_max); % what node #
                ksym = sympmap(i,j,j_max);
                k_e = k + 1;
                k_w = k - 1;
                k_s = k + j_max;
                k_n = ksym + 1;

                thisk_k = M{mat(k)}.k;
                thisk_ke = M{mat(k_e)}.k;
                thisk_kw = M{mat(k_w)}.k;
                thisk_kn = M{mat(k_n)}.k;
                thisk_ks = M{mat(k_s)}.k;
                % thisrhoc_p = M{mat(k)}.rhoCp;

                if any(MODindeces == mat(k)) % If coolant
                    heatremoval = MOD_rhoc_p*dTdz*w(k);
                    if mat(k_e) == 1
                        flux_e = h/Deltax;
                    else
                        flux_e = (thisk_ke+thisk_k)/Deltax^2;
                    end
                    if mat(k_w) == 1
                        flux_w = h/Deltax;
                    else
                        flux_w = (thisk_kw+thisk_k)/Deltax^2;
                    end
                    if mat(k_n) == 1
                        flux_n = h/Deltay;
                    else
                        flux_n = (thisk_kn+thisk_k)/Deltay^2;
                    end
                    if mat(k_s) == 1
                        flux_s = h/Deltay;
                    else
                        flux_s = (thisk_ks+thisk_k)/Deltay^2;
                    end
                else % Fuel
                    heatremoval = modCorr*MOD_rhoc_p*dTdz*w(k);
                    if any(MODindeces == mat(k_e))
                        flux_e = h/Deltax;
                    else
                        flux_e = (thisk_ke+thisk_k)/Deltax^2;
                    end
                    if any(MODindeces == mat(k_w))
                        flux_w = h/Deltax;
                    else
                        flux_w = (thisk_kw+thisk_k)/Deltax^2;
                    end
                    if any(MODindeces == mat(k_n))
                        flux_n = h/Deltay;
                    else
                        flux_n = (thisk_kn+thisk_k)/Deltay^2;
                    end
                    if any(MODindeces == mat(k_s))
                        flux_s = h/Deltay;
                    else
                        flux_s = (thisk_ks+thisk_k)/Deltay^2;
                    end
                end
                A(k,k) = flux_e + flux_w + flux_n + flux_s;
                A(k,k_e) = A(k,k_e) - flux_e;
                A(k,k_w) = A(k,k_w) - flux_w;
                A(k,k_n) = A(k,k_n) - flux_n;
                A(k,k_s) = A(k,k_s) - flux_s;
                Q(k) = (q3prime(k) - heatremoval);
            end
        end

        % Center Boundary
        i = 1;
        j = 1;
        k = pmap(i, j, j_max); % what node #
        k_e = k + 1;
        k_s = k + i_max;
        k_n = k_e;
        k_w = k_s;

        thisk_k = M{mat(k)}.k;
        thisk_ke = M{mat(k_e)}.k;
        thisk_kw = M{mat(k_w)}.k;
        thisk_kn = M{mat(k_n)}.k;
        thisk_ks = M{mat(k_s)}.k;
        % thisrhoc_p = M{mat(k)}.rhoCp;

        if any(MODindeces == mat(k)) % If coolant
            heatremoval = MOD_rhoc_p*dTdz*w(k);
            if mat(k_e) == 1
                flux_e = h/Deltax;
            else
                flux_e = (thisk_ke+thisk_k)/Deltax^2;
            end
            if mat(k_w) == 1
                flux_w = h/Deltax;
            else
                flux_w = (thisk_kw+thisk_k)/Deltax^2;
            end
            if mat(k_n) == 1
                flux_n = h/Deltay;
            else
                flux_n = (thisk_kn+thisk_k)/Deltay^2;
            end
            if mat(k_s) == 1
                flux_s = h/Deltay;
            else
                flux_s = (thisk_ks+thisk_k)/Deltay^2;
            end
        else % Fuel
            heatremoval = modCorr*MOD_rhoc_p*dTdz*w(k);
            if any(MODindeces == mat(k_e))
                flux_e = h/Deltax;
            else
                flux_e = (thisk_ke+thisk_k)/Deltax^2;
            end
            if any(MODindeces == mat(k_w))
                flux_w = h/Deltax;
            else
                flux_w = (thisk_kw+thisk_k)/Deltax^2;
            end
            if any(MODindeces == mat(k_n))
                flux_n = h/Deltay;
            else
                flux_n = (thisk_kn+thisk_k)/Deltay^2;
            end
            if any(MODindeces == mat(k_s))
                flux_s = h/Deltay;
            else
                flux_s = (thisk_ks+thisk_k)/Deltay^2;
            end
        end
        A(k,k) = flux_e + flux_w + flux_n + flux_s;
        A(k,k_e) = A(k,k_e) - flux_e;
        A(k,k_w) = A(k,k_w) - flux_w;
        A(k,k_n) = A(k,k_n) - flux_n;
        A(k,k_s) = A(k,k_s) - flux_s;
        Q(k) = (q3prime(k) - heatremoval);

        fprintf('Solving for %i degrees of freedom...\n', i_max*j_max);
        T = A\Q; % Axial Flow Field [cm/s] -> forms a heat removal term
        clearvars A Q
        fprintf('Complete!\n');
        TPlot = reshape(T, i_max, j_max);
        figure(1);
        surf(x, y, TPlot)
        ylabel('y');
        xlabel('x');
        title('Temperature Surface');
        drawnow;
    end

%% Iterative Solver: Check Solution + Update Boron
% ------------------------------------------------------------------------------

    keff_iter = keff;

    if keff_iter < 1.0
        myBORLo = myBOR;
    elseif keff_iter > 1.0
        myBORHi = myBOR;
    end

end

%% Final Plot
% ------------------------------------------------------------------------------

Tplot = reshape(T, i_max, j_max);
Tref1 = rot90(Tplot,3);
Tref1 = Tref1(:,1:end-1);
Tref1Plot = horzcat(Tref1, Tplot);
Tref2 = rot90(Tref1Plot, 2);
Tref2 = Tref2(1:end-1,:);
Tref2Plot = vertcat(Tref2, Tref1Plot);

figure(1);
% Plot temp surface
surf(fullx,fully,Tref2Plot);
ylabel('y [cm]');
xlabel('x [cm]');
zlabel('Temperature [K]')
title('Steady-State Solution for a Reactor');
saveas(figure(1),fullfile(myCWD,subfolder,'tempSurface.jpg'));

%% Save Results
% ------------------------------------------------------------------------------

resultmat = 'myVals.mat';
resultOut = fullfile(myCWD,subfolder,resultmat);

% Parameters of Interest
% myBOR [ppm]
% keff_iter []
save(resultOut, 'TM_nom', 'myBOR', 'keff_iter');

% And Solution Matrices
save(fullfile(myCWD,subfolder,'w.mat'), 'w');
save(fullfile(myCWD,subfolder,'phi.mat'), 'phi');
save(fullfile(myCWD,subfolder,'xPlot.mat'), 'fullx');
save(fullfile(myCWD,subfolder,'yPlot.mat'), 'fully');
save(fullfile(myCWD,subfolder,'TPlot.mat'), 'Tref2Plot');

%% Functions
% ------------------------------------------------------------------------------



function D = myD(tfu, tmo, bor, g, matno)
    global TFs;
    global TMs;
    global BORs;
    global D1;
    global D2;
    if g==1
        D = interpn(TFs{matno}, TMs{matno}, BORs{matno}, D1{matno}, tfu, tmo, bor,'spline');
    elseif g==2
        D = interpn(TFs{matno}, TMs{matno}, BORs{matno}, D2{matno}, tfu, tmo, bor,'spline');
    else
        D = 0;
    end
end

function Sigma_R = mySigma_R(tfu, tmo, bor, g, matno)
    global TFs;
    global TMs;
    global BORs;
    global Sigma_R1;
    global Sigma_R2;
    if g==1
        Sigma_R = interpn(TFs{matno}, TMs{matno}, BORs{matno}, Sigma_R1{matno}, tfu, tmo, bor,'spline');
    elseif g==2
        Sigma_R = interpn(TFs{matno}, TMs{matno}, BORs{matno}, Sigma_R2{matno}, tfu, tmo, bor,'spline');
    else
        Sigma_R = 0;
    end
end

function Sigma_s =mySigma_s(tfu, tmo, bor, gprime, g, matno)
    global TFs;
    global TMs;
    global BORs;
    % global Sigma_s11;
    global Sigma_s12;
    global Sigma_s21;
    % global Sigma_s22;
    if gprime==1
        if g==1
            Sigma_s = 0; % interpn(FUL_TFs, FUL_TMs, FUL_BORs, FUL_Sigma_s11, tfu, tmo, bor);
        elseif g==2
            Sigma_s = interpn(TFs{matno}, TMs{matno}, BORs{matno}, Sigma_s12{matno}, tfu, tmo, bor,'spline');
        end
    elseif gprime==2
        if g==1
            Sigma_s = interpn(TFs{matno}, TMs{matno}, BORs{matno}, Sigma_s21{matno}, tfu, tmo, bor,'spline');
        elseif g==2
            Sigma_s = 0; % interpn(FUL_TFs, FUL_TMs, FUL_BORs, FUL_Sigma_s22, tfu, tmo, bor);
        end
    else
        Sigma_s = 0;
    end
end

function nuSigma_f = mynuSigma_f(tfu, tmo, bor, g, matno)
    global TFs;
    global TMs;
    global BORs;
    global nuSigma_f1;
    global nuSigma_f2;
    if g==1
        nuSigma_f = interpn(TFs{matno}, TMs{matno}, BORs{matno}, nuSigma_f1{matno}, tfu, tmo, bor,'spline');
    elseif g==2
        nuSigma_f = interpn(TFs{matno}, TMs{matno}, BORs{matno}, nuSigma_f2{matno}, tfu, tmo, bor,'spline');
    else
        nuSigma_f = 0;
    end
end

% Pointer mapping function for 2D->1D transform of solution
function k = pmap(i, j, j_max)
    k = j + (i-1)*j_max;
end

function symk = sympmap(i, j, j_max)
    symk = i + (j-1)*j_max;
end