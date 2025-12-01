%% 2D Coupled SMR
% ------------------------------------------------------------------------------
% Adam Buchalter
%
% desc
% ------------------------------------------------------------------------------
clear all; close all; clc;

%% Import Data
% ------------------------------------------------------------------------------
myCWD = pwd;

% Thermal Data
FUL = readtable(fullfile(myCWD,'proj6Data\\thermalData\\Fuel_Lump.csv'));
MOD = readtable(fullfile(myCWD,'proj6Data\\thermalData\\H2O_Lump.csv'));
% Cell Array of Materials
M = {FUL, MOD};

% Neutron Data
myMat = load(fullfile(myCWD,'proj6Data\\neutronData\\interpTables.mat'));
global FUL_TFs;
global FUL_TMs;
global FUL_BORs;
global FUL_D1;
global FUL_D2;
global FUL_Sigma_R1;
global FUL_Sigma_R2;
global FUL_Sigma_a1;
global FUL_Sigma_a2;
global FUL_Sigma_s11;
global FUL_Sigma_s12;
global FUL_Sigma_s21;
global FUL_Sigma_s22;
global FUL_nuSigma_f1;
global FUL_nuSigma_f2;
FUL_TFs = myMat.TFs;
FUL_TMs = myMat.TMs;
FUL_BORs = myMat.BORs;
FUL_D1 = myMat.D1;
FUL_D2 = myMat.D2;
FUL_Sigma_R1 = myMat.Sigma_R1;
FUL_Sigma_R2 = myMat.Sigma_R2;
FUL_Sigma_a1 = myMat.Sigma_a1;
FUL_Sigma_a2 = myMat.Sigma_a2;
FUL_Sigma_s11 = myMat.Sigma_s11;
FUL_Sigma_s12 = myMat.Sigma_s12;
FUL_Sigma_s21 = myMat.Sigma_s21;
FUL_Sigma_s22 = myMat.Sigma_s22;
FUL_nuSigma_f1 = myMat.nuSigma_f1;
FUL_nuSigma_f2 = myMat.nuSigma_f2;
CHI = [1 0];

myMat = load(fullfile(myCWD,'proj6Data\\neutronData\\interpTablesMod.mat'));
global MOD_TFs;
global MOD_TMs;
global MOD_BORs;
global MOD_D1;
global MOD_D2;
global MOD_Sigma_R1;
global MOD_Sigma_R2;
global MOD_Sigma_a1;
global MOD_Sigma_a2;
global MOD_Sigma_s11;
global MOD_Sigma_s12;
global MOD_Sigma_s21;
global MOD_Sigma_s22;
global MOD_nuSigma_f1;
global MOD_nuSigma_f2;
MOD_TFs = myMat.TFs;
MOD_TMs = myMat.TMs;
MOD_BORs = myMat.BORs;
MOD_D1 = myMat.D1;
MOD_D2 = myMat.D2;
MOD_Sigma_R1 = myMat.Sigma_R1;
MOD_Sigma_R2 = myMat.Sigma_R2;
MOD_Sigma_a1 = myMat.Sigma_a1;
MOD_Sigma_a2 = myMat.Sigma_a2;
MOD_Sigma_s11 = myMat.Sigma_s11;
MOD_Sigma_s12 = myMat.Sigma_s12;
MOD_Sigma_s21 = myMat.Sigma_s21;
MOD_Sigma_s22 = myMat.Sigma_s22;
MOD_nuSigma_f1 = myMat.nuSigma_f1;
MOD_nuSigma_f2 = myMat.nuSigma_f2;

% Core Layout
% LAYOUT = [
%     1 1
%     1 1
% ];
LAYOUT = [
    1 1 1 1 2
    1 1 1 1 2
    1 1 1 2 2
    1 1 2 2 2
    2 2 2 2 2
];

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


TF_nom = 850.0; % Nominal Fuel Temperature [K]
TM_nom = 557.0; % Nominal Moderator Temperature [K]
TFUL = TF_nom;
TMOD = TM_nom;
myBOR = 600.0; % Init Soluble Poison Concentration [ppm]

%% Intelligently Find Smallest Diffusion Length for Resolution-Setting
% ------------------------------------------------------------------------------
% The smallest diffusion length is important because if our cells are larger
% than the diffusion length for a particular group then we cannot reliably
% capture transport phenomena

L = 50; % [cm] Arbitrarily Large Initial Guess
for i = 1:size(FUL_TFs)
    for j = 1:size(FUL_TMs)
        for k = 1:size(FUL_BORs)
            thisL = sqrt(FUL_D1(i,j,k)/FUL_Sigma_a1(i,j,k)); % Characteristic Diffusion Length [cm]
            if thisL < L % If smaller, use as new normalization length
                L = thisL;
            end
            thisL = sqrt(FUL_D2(i,j,k)/FUL_Sigma_a2(i,j,k)); % Characteristic Diffusion Length [cm]
            if thisL < L % If smaller, use as new normalization length
                L = thisL;
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

% Define variables for power-iteration
residual = 1.0E5; % init residual
epsilon = 1.0E-16; % drive residual down to this value before terminating

% % Define time stepping
% Deltat = 1; % [s]

%% Nondimensional Domain Prep
% ------------------------------------------------------------------------------

% Define physical domain
size = Deltax*(i_max-1)*2; % FULL Domain Size [cm]

% Unitless Constants
T_r = T_infty; % [K]

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

% q3prime = zeros(i_max*j_max,1);

% File Info
subfolder='results\\project6\\';
mkdir(fullfile(myCWD,subfolder));

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
        if thisDomain == 1 % Fuel
            A(k, k) = 1.0;
            S(k) = myVel; % Bulk Core Velocity Source [cm/s]
        end
        if thisDomain == 2 % Moderator
            A(k,k) = -2.0*(factorx + factory);
            A(k,k_e) = factorx;
            A(k,k_w) = factorx;
            A(k,k_n) = factory;
            A(k,k_s) = factory;
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
        if thisDomain == 1 % Fuel
            A(k, k) = 1.0;
            S(k) = myVel; % Bulk Core Velocity Source [cm/s]
        end
        if thisDomain == 2 % Moderator
            A(k,k) = -2.0*(factorx + factory);
            A(k,k_e) = A(k,k_e) + factorx;
            A(k,k_w) = A(k,k_w) + factorx;
            A(k,k_n) = A(k,k_n) + factory;
            A(k,k_s) = A(k,k_s) + factory;
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
        if thisDomain == 1 % Fuel
            A(k, k) = 1.0;
            S(k) = myVel; % Bulk Core Velocity Source [cm/s]
        end
        if thisDomain == 2 % Moderator
            A(k,k) = -2.0*(factorx + factory);
            A(k,k_e) = A(k,k_e) + factorx;
            A(k,k_w) = A(k,k_w) + factorx;
            A(k,k_n) = A(k,k_n) + factory;
            A(k,k_s) = A(k,k_s) + factory;
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
if thisDomain == 1 % Fuel
    A(k, k) = 1.0;
    S(k) = myVel; % Bulk Core Velocity Source [cm/s]
end
if thisDomain == 2 % Moderator
    A(k,k) = -2.0*(factorx + factory);
    A(k,k_e) = A(k,k_e) + factorx;
    A(k,k_w) = A(k,k_w) + factorx;
    A(k,k_n) = A(k,k_n) + factory;
    A(k,k_s) = A(k,k_s) + factory;
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

%% Iterative Solver: Neutronics
% ------------------------------------------------------------------------------
fprintf('\n===ITERATIVE SOLVER===\n');
fprintf('Initializing...\n');
% Init Thermal Matrices
A = spalloc(i_max*j_max, i_max*j_max, 5*i_max*j_max); % Sparsely allocate Diffusion Operator with 5 bands
Q = zeros(i_max*j_max, 1); % Allocate Heat Source/Removal Vector
% Init Neutronics Matrices
H = spalloc(G*i_max*j_max, G*i_max*j_max, 5*G*i_max*j_max); % Sparsely allocate Streaming/Absorption Operator with 5 bands for each energy group
S = spalloc(G*i_max*j_max, G*i_max*j_max, (G-1)*i_max*j_max); % Sparsely allocate Scattering Source Operator with bands for up/downscatter
F = spalloc(G*i_max*j_max, G*i_max*j_max, 1*G*i_max*j_max); % Sparsely allocate Fission Source Operator with 1 band for each energy group

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

fprintf('\n---Solving for Neutron Flux Profile---\n');

% Define Matrices
fprintf('Building Matrices...\n');
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

            if thisMat == 1 % Fuel
                thisD = FUL_D(thisT, TMOD, myBOR, group);
                thisSig_R = FUL_Sigma_R(thisT, TMOD, myBOR, group);
            elseif thisMat == 2 % Moderator
                thisD = MOD_D(TFUL, thisT, myBOR, group);
                thisSig_R = MOD_Sigma_R(TFUL, thisT, myBOR, group);
            end

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
                
                if thisMat == 1 % Fuel
                    nuSig_fgprime = FUL_nuSigma_f(thisT, TMOD, myBOR, gprime);
                    Sig_sgprimeg = FUL_Sigma_s(thisT, TMOD, myBOR, gprime, group);
                elseif thisMat == 2 % Moderator
                    nuSig_fgprime = MOD_nuSigma_f(TFUL, thisT, myBOR, gprime);
                    Sig_sgprimeg = MOD_Sigma_s(TFUL, thisT, myBOR, gprime, group);
                end 

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

            if thisMat == 1 % Fuel
                thisD = FUL_D(thisT, TMOD, myBOR, group);
                thisSig_R = FUL_Sigma_R(thisT, TMOD, myBOR, group);
            elseif thisMat == 2 % Moderator
                thisD = MOD_D(TFUL, thisT, myBOR, group);
                thisSig_R = MOD_Sigma_R(TFUL, thisT, myBOR, group);
            end

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
                
                if thisMat == 1 % Fuel
                    nuSig_fgprime = FUL_nuSigma_f(thisT, TMOD, myBOR, gprime);
                    Sig_sgprimeg = FUL_Sigma_s(thisT, TMOD, myBOR, gprime, group);
                elseif thisMat == 2 % Moderator
                    nuSig_fgprime = MOD_nuSigma_f(TFUL, thisT, myBOR, gprime);
                    Sig_sgprimeg = MOD_Sigma_s(TFUL, thisT, myBOR, gprime, group);
                end 

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

            if thisMat == 1 % Fuel
                thisD = FUL_D(thisT, TMOD, myBOR, group);
                thisSig_R = FUL_Sigma_R(thisT, TMOD, myBOR, group);
            elseif thisMat == 2 % Moderator
                thisD = MOD_D(TFUL, thisT, myBOR, group);
                thisSig_R = MOD_Sigma_R(TFUL, thisT, myBOR, group);
            end

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
                
                if thisMat == 1 % Fuel
                    nuSig_fgprime = FUL_nuSigma_f(thisT, TMOD, myBOR, gprime);
                    Sig_sgprimeg = FUL_Sigma_s(thisT, TMOD, myBOR, gprime, group);
                elseif thisMat == 2 % Moderator
                    nuSig_fgprime = MOD_nuSigma_f(TFUL, thisT, myBOR, gprime);
                    Sig_sgprimeg = MOD_Sigma_s(TFUL, thisT, myBOR, gprime, group);
                end 

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

    if thisMat == 1 % Fuel
        thisD = FUL_D(thisT, TMOD, myBOR, group);
        thisSig_R = FUL_Sigma_R(thisT, TMOD, myBOR, group);
    elseif thisMat == 2 % Moderator
        thisD = MOD_D(TFUL, thisT, myBOR, group);
        thisSig_R = MOD_Sigma_R(TFUL, thisT, myBOR, group);
    end

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
        
        if thisMat == 1 % Fuel
            nuSig_fgprime = FUL_nuSigma_f(thisT, TMOD, myBOR, gprime);
            Sig_sgprimeg = FUL_Sigma_s(thisT, TMOD, myBOR, gprime, group);
        elseif thisMat == 2 % Moderator
            nuSig_fgprime = MOD_nuSigma_f(TFUL, thisT, myBOR, gprime);
            Sig_sgprimeg = MOD_Sigma_s(TFUL, thisT, myBOR, gprime, group);
        end 

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

q3prime = phi(1:i_max*j_max)+phi(1+i_max*j_max:2*i_max*j_max);
q3prime = (totLinPwr/4)*q3prime/(sum(q3prime,'all')*Deltax*Deltay);
% q3prime = (totPwr/4)*q3prime/(sum(q3prime,'all')*Deltax*Deltay);

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
% Q = q3prime - MOD_rhoc_p*dTdz*w;
for i = 2:i_max-1
    for j = 2:j_max-1
        k = pmap(i, j, j_max); % what node #
        k_e = k + 1;
        k_w = k - 1;
        k_n = k - j_max;
        k_s = k + j_max;

        % thisMat = Domain(k); % what material
        % matnum = mat(k); % what material
        thisk_k = M{mat(k)}.k;
        thisk_ke = M{mat(k_e)}.k;
        thisk_kw = M{mat(k_w)}.k;
        thisk_kn = M{mat(k_n)}.k;
        thisk_ks = M{mat(k_s)}.k;
        % thisrhoc_p = M{mat(k)}.rhoCp;

        if mat(k) == 1 % Fuel
            heatremoval = modCorr*MOD_rhoc_p*dTdz*w(k);
            if mat(k_e) == 2
                flux_e = h/Deltax;
                % heatremoval = heatremoval - h*T_infty/Deltax;
            else
                flux_e = (thisk_ke+thisk_k)/Deltax^2;
            end
            if mat(k_w) == 2
                flux_w = h/Deltax;
                % heatremoval = heatremoval - h*T_infty/Deltax;
            else
                flux_w = (thisk_kw+thisk_k)/Deltax^2;
            end
            if mat(k_n) == 2
                flux_n = h/Deltay;
                % heatremoval = heatremoval - h*T_infty/Deltay;
            else
                flux_n = (thisk_kn+thisk_k)/Deltay^2;
            end
            if mat(k_s) == 2
                flux_s = h/Deltay;
                % heatremoval = heatremoval - h*T_infty/Deltay;
            else
                flux_s = (thisk_ks+thisk_k)/Deltay^2;
            end
        else
            % A(k,k) = 1.0;
            % Q(k) = T_infty;
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
        end
        A(k,k) = flux_e + flux_w + flux_n + flux_s;
        A(k,k_e) = A(k,k_e) - flux_e;
        A(k,k_w) = A(k,k_w) - flux_w;
        A(k,k_n) = A(k,k_n) - flux_n;
        A(k,k_s) = A(k,k_s) - flux_s;
        Q(k) = (q3prime(k) - heatremoval);

        % pointer mapping goes row-by-row to assemble Coeff. Matrix
        % A(k,k) = (thisk_ke+2.0*thisk_k+thisk_kw)/Deltax^2 + (thisk_kn+2.0*thisk_k+thisk_ks)/Deltay^2;
        % A(k,k_e) = -(thisk_ke+thisk_k)/Deltax^2;
        % A(k,k_w) = -(thisk_kw+thisk_k)/Deltax^2;
        % A(k,k_n) = -(thisk_kn+thisk_k)/Deltay^2;
        % A(k,k_s) = -(thisk_ks+thisk_k)/Deltay^2;
        % %
        % Q(k) = (q3prime(k) - heatremoval)/1E4;
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

        if mat(k) == 1 % Fuel
            heatremoval = modCorr*MOD_rhoc_p*dTdz*w(k);
            if mat(k_e) == 2
                flux_e = h/Deltax;
                % heatremoval = heatremoval - h*T_infty/Deltax;
            else
                flux_e = (thisk_ke+thisk_k)/Deltax^2;
            end
            if mat(k_w) == 2
                flux_w = h/Deltax;
                % heatremoval = heatremoval - h*T_infty/Deltax;
            else
                flux_w = (thisk_kw+thisk_k)/Deltax^2;
            end
            if mat(k_n) == 2
                flux_n = h/Deltay;
                % heatremoval = heatremoval - h*T_infty/Deltay;
            else
                flux_n = (thisk_kn+thisk_k)/Deltay^2;
            end
            if mat(k_s) == 2
                flux_s = h/Deltay;
                % heatremoval = heatremoval - h*T_infty/Deltay;
            else
                flux_s = (thisk_ks+thisk_k)/Deltay^2;
            end
        else
            % A(k,k) = 1.0;
            % Q(k) = T_infty;
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
        end
        A(k,k) = flux_e + flux_w + flux_n + flux_s;
        A(k,k_e) = A(k,k_e) - flux_e;
        A(k,k_w) = A(k,k_w) - flux_w;
        A(k,k_n) = A(k,k_n) - flux_n;
        A(k,k_s) = A(k,k_s) - flux_s;
        Q(k) = (q3prime(k) - heatremoval);

        % pointer mapping goes row-by-row to assemble Coeff. Matrix
        % A(k,k) = (thisk_ke+2.0*thisk_k+thisk_kw)/Deltax^2 + (thisk_kn+2.0*thisk_k+thisk_ks)/Deltay^2;
        % A(k,k_e) = A(k,k_e) - (thisk_ke+thisk_k)/Deltax^2;
        % A(k,k_w) = A(k,k_w) - (thisk_kw+thisk_k)/Deltax^2;
        % A(k,k_n) = A(k,k_n) - (thisk_kn+thisk_k)/Deltay^2;
        % A(k,k_s) = A(k,k_s) - (thisk_ks+thisk_k)/Deltay^2;
        %
        % Q(k) = (q3prime(k) - heatremoval)/1E4;
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

        % thisMat = Domain(k); % what material
        % matnum = mat(k); % what material
        thisk_k = M{mat(k)}.k;
        thisk_ke = M{mat(k_e)}.k;
        thisk_kw = M{mat(k_w)}.k;
        thisk_kn = M{mat(k_n)}.k;
        thisk_ks = M{mat(k_s)}.k;
        % thisrhoc_p = M{mat(k)}.rhoCp;

        if mat(k) == 1 % Fuel
            heatremoval = modCorr*MOD_rhoc_p*dTdz*w(k);
            if mat(k_e) == 2
                flux_e = h/Deltax;
                % heatremoval = heatremoval - h*T_infty/Deltax;
            else
                flux_e = (thisk_ke+thisk_k)/Deltax^2;
            end
            if mat(k_w) == 2
                flux_w = h/Deltax;
                % heatremoval = heatremoval - h*T_infty/Deltax;
            else
                flux_w = (thisk_kw+thisk_k)/Deltax^2;
            end
            if mat(k_n) == 2
                flux_n = h/Deltay;
                % heatremoval = heatremoval - h*T_infty/Deltay;
            else
                flux_n = (thisk_kn+thisk_k)/Deltay^2;
            end
            if mat(k_s) == 2
                flux_s = h/Deltay;
                % heatremoval = heatremoval - h*T_infty/Deltay;
            else
                flux_s = (thisk_ks+thisk_k)/Deltay^2;
            end
        else
            % A(k,k) = 1.0;
            % Q(k) = T_infty;
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
        end
        A(k,k) = flux_e + flux_w + flux_n + flux_s;
        A(k,k_e) = A(k,k_e) - flux_e;
        A(k,k_w) = A(k,k_w) - flux_w;
        A(k,k_n) = A(k,k_n) - flux_n;
        A(k,k_s) = A(k,k_s) - flux_s;
        Q(k) = (q3prime(k) - heatremoval);

        % pointer mapping goes row-by-row to assemble Coeff. Matrix
        % A(k,k) = (thisk_ke+2.0*thisk_k+thisk_kw)/Deltax^2 + (thisk_kn+2.0*thisk_k+thisk_ks)/Deltay^2;
        % A(k,k_e) = A(k,k_e) - (thisk_ke+thisk_k)/Deltax^2;
        % A(k,k_w) = A(k,k_w) - (thisk_kw+thisk_k)/Deltax^2;
        % A(k,k_n) = A(k,k_n) - (thisk_kn+thisk_k)/Deltay^2;
        % A(k,k_s) = A(k,k_s) - (thisk_ks+thisk_k)/Deltay^2;
        %
        % Q(k) = (q3prime(k) - heatremoval)/1E4;
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

% thisMat = Domain(k); % what material
% matnum = mat(k); % what material
thisk_k = M{mat(k)}.k;
thisk_ke = M{mat(k_e)}.k;
thisk_kw = M{mat(k_w)}.k;
thisk_kn = M{mat(k_n)}.k;
thisk_ks = M{mat(k_s)}.k;
% thisrhoc_p = M{mat(k)}.rhoCp;

if mat(k) == 1 % Fuel
    heatremoval = modCorr*MOD_rhoc_p*dTdz*w(k);
    if mat(k_e) == 2
        flux_e = h/Deltax;
        % heatremoval = heatremoval - h*T_infty/Deltax;
    else
        flux_e = (thisk_ke+thisk_k)/Deltax^2;
    end
    if mat(k_w) == 2
        flux_w = h/Deltax;
        % heatremoval = heatremoval - h*T_infty/Deltax;
    else
        flux_w = (thisk_kw+thisk_k)/Deltax^2;
    end
    if mat(k_n) == 2
        flux_n = h/Deltay;
        % heatremoval = heatremoval - h*T_infty/Deltay;
    else
        flux_n = (thisk_kn+thisk_k)/Deltay^2;
    end
    if mat(k_s) == 2
        flux_s = h/Deltay;
        % heatremoval = heatremoval - h*T_infty/Deltay;
    else
        flux_s = (thisk_ks+thisk_k)/Deltay^2;
    end
else
    % A(k,k) = 1.0;
    % Q(k) = T_infty;
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
end
A(k,k) = flux_e + flux_w + flux_n + flux_s;
A(k,k_e) = A(k,k_e) - flux_e;
A(k,k_w) = A(k,k_w) - flux_w;
A(k,k_n) = A(k,k_n) - flux_n;
A(k,k_s) = A(k,k_s) - flux_s;
Q(k) = (q3prime(k) - heatremoval);

% pointer mapping goes row-by-row to assemble Coeff. Matrix
% A(k,k) = (thisk_ke+2.0*thisk_k+thisk_kw)/Deltax^2 + (thisk_kn+2.0*thisk_k+thisk_ks)/Deltay^2;
% A(k,k_e) = A(k,k_e) - (thisk_ke+thisk_k)/Deltax^2;
% A(k,k_w) = A(k,k_w) - (thisk_kw+thisk_k)/Deltax^2;
% A(k,k_n) = A(k,k_n) - (thisk_kn+thisk_k)/Deltay^2;
% A(k,k_s) = A(k,k_s) - (thisk_ks+thisk_k)/Deltay^2;
% %
% Q(k) = (q3prime(k) - heatremoval)/1E4;

fprintf('Solving for %i degrees of freedom...\n', i_max*j_max);
T = A\Q; % Axial Flow Field [cm/s] -> forms a heat removal term
fprintf('Complete!\n');
TPlot = reshape(T, i_max, j_max);
figure(1);
surf(x, y, TPlot)
ylabel('y');
xlabel('x');
title('Temperature Surface');
drawnow;

%% Functions
% ------------------------------------------------------------------------------

function D = FUL_D(tfu, tmo, bor, g)
    global FUL_TFs;
    global FUL_TMs;
    global FUL_BORs;
    global FUL_D1;
    global FUL_D2;
    if g==1
        D = interpn(FUL_TFs, FUL_TMs, FUL_BORs, FUL_D1, tfu, tmo, bor);
    elseif g==2
        D = interpn(FUL_TFs, FUL_TMs, FUL_BORs, FUL_D2, tfu, tmo, bor);
    else
        D = 0;
    end
end

function Sigma_R = FUL_Sigma_R(tfu, tmo, bor, g)
    global FUL_TFs;
    global FUL_TMs;
    global FUL_BORs;
    global FUL_Sigma_R1;
    global FUL_Sigma_R2;
    if g==1
        Sigma_R = interpn(FUL_TFs, FUL_TMs, FUL_BORs, FUL_Sigma_R1, tfu, tmo, bor);
    elseif g==2
        Sigma_R = interpn(FUL_TFs, FUL_TMs, FUL_BORs, FUL_Sigma_R2, tfu, tmo, bor);
    else
        Sigma_R = 0;
    end
end

% function Sigma_a = FUL_Sigma_a(tfu, tmo, bor, g)
%     if g==1
%         Sigma_a = interpn(FUL_TFs, FUL_TMs, FUL_BORs, FUL_Sigma_a1, tfu, tmo, bor);
%     elseif g==2
%         Sigma_a = interpn(FUL_TFs, FUL_TMs, FUL_BORs, FUL_Sigma_a2, tfu, tmo, bor);
%     else
%         Sigma_a = 0;
%     end
% end

function Sigma_s = FUL_Sigma_s(tfu, tmo, bor, gprime, g)
    global FUL_TFs;
    global FUL_TMs;
    global FUL_BORs;
    global FUL_Sigma_s11;
    global FUL_Sigma_s12;
    global FUL_Sigma_s21;
    global FUL_Sigma_s22;
    if gprime==1
        if g==1
            Sigma_s = 0; % interpn(FUL_TFs, FUL_TMs, FUL_BORs, FUL_Sigma_s11, tfu, tmo, bor);
        elseif g==2
            Sigma_s = interpn(FUL_TFs, FUL_TMs, FUL_BORs, FUL_Sigma_s12, tfu, tmo, bor);
        end
    elseif gprime==2
        if g==1
            Sigma_s = interpn(FUL_TFs, FUL_TMs, FUL_BORs, FUL_Sigma_s21, tfu, tmo, bor);
        elseif g==2
            Sigma_s = 0; % interpn(FUL_TFs, FUL_TMs, FUL_BORs, FUL_Sigma_s22, tfu, tmo, bor);
        end
    else
        Sigma_s = 0;
    end
end

function nuSigma_f = FUL_nuSigma_f(tfu, tmo, bor, g)
    global FUL_TFs;
    global FUL_TMs;
    global FUL_BORs;
    global FUL_nuSigma_f1;
    global FUL_nuSigma_f2;
    if g==1
        nuSigma_f = interpn(FUL_TFs, FUL_TMs, FUL_BORs, FUL_nuSigma_f1, tfu, tmo, bor);
    elseif g==2
        nuSigma_f = interpn(FUL_TFs, FUL_TMs, FUL_BORs, FUL_nuSigma_f2, tfu, tmo, bor);
    else
        nuSigma_f = 0;
    end
end

function D = MOD_D(tfu, tmo, bor, g)
    global MOD_TFs;
    global MOD_TMs;
    global MOD_BORs;
    global MOD_D1;
    global MOD_D2;
    if g==1
        D = interpn(MOD_TFs, MOD_TMs, MOD_BORs, MOD_D1, tfu, tmo, bor);
    elseif g==2
        D = interpn(MOD_TFs, MOD_TMs, MOD_BORs, MOD_D2, tfu, tmo, bor);
    else
        D = 0;
    end
end

function Sigma_R = MOD_Sigma_R(tfu, tmo, bor, g)
    global MOD_TFs;
    global MOD_TMs;
    global MOD_BORs;
    global MOD_Sigma_R1;
    global MOD_Sigma_R2;
    if g==1
        Sigma_R = interpn(MOD_TFs, MOD_TMs, MOD_BORs, MOD_Sigma_R1, tfu, tmo, bor);
    elseif g==2
        Sigma_R = interpn(MOD_TFs, MOD_TMs, MOD_BORs, MOD_Sigma_R2, tfu, tmo, bor);
    else
        Sigma_R = 0;
    end
end

% function Sigma_a = MOD_Sigma_a(tfu, tmo, bor, g)
%     if g==1
%         Sigma_a = interpn(MOD_TFs, MOD_TMs, MOD_BORs, MOD_Sigma_a1, tfu, tmo, bor);
%     elseif g==2
%         Sigma_a = interpn(MOD_TFs, MOD_TMs, MOD_BORs, MOD_Sigma_a2, tfu, tmo, bor);
%     else
%         Sigma_a = 0;
%     end
% end

function Sigma_s = MOD_Sigma_s(tfu, tmo, bor, gprime, g)
    global MOD_TFs;
    global MOD_TMs;
    global MOD_BORs;
    global MOD_Sigma_s11;
    global MOD_Sigma_s12;
    global MOD_Sigma_s21;
    global MOD_Sigma_s22;
    if gprime==1
        if g==1
            Sigma_s = 0; % interpn(MOD_TFs, MOD_TMs, MOD_BORs, MOD_Sigma_s11, tfu, tmo, bor);
        elseif g==2
            Sigma_s = interpn(MOD_TFs, MOD_TMs, MOD_BORs, MOD_Sigma_s12, tfu, tmo, bor);
        end
    elseif gprime==2
        if g==1
            Sigma_s = interpn(MOD_TFs, MOD_TMs, MOD_BORs, MOD_Sigma_s21, tfu, tmo, bor);
        elseif g==2
            Sigma_s = 0; %  interpn(MOD_TFs, MOD_TMs, MOD_BORs, MOD_Sigma_s22, tfu, tmo, bor);
        end
    else
        Sigma_s = 0;
    end
end

function nuSigma_f = MOD_nuSigma_f(tfu, tmo, bor, g)
    global MOD_TFs;
    global MOD_TMs;
    global MOD_BORs;
    global MOD_nuSigma_f1;
    global MOD_nuSigma_f2;
    if g==1
        nuSigma_f = interpn(MOD_TFs, MOD_TMs, MOD_BORs, MOD_nuSigma_f1, tfu, tmo, bor);
    elseif g==2
        nuSigma_f = interpn(MOD_TFs, MOD_TMs, MOD_BORs, MOD_nuSigma_f2, tfu, tmo, bor);
    else
        nuSigma_f = 0;
    end
end

% Pointer mapping function for 2D->1D transform of solution
function k = pmap(i, j, j_max)
    k = j + (i-1)*j_max;
    % k = i + (j-1)*i_max;
end

function symk = sympmap(i, j, j_max)
    symk = i + (j-1)*j_max;
    % symk = (j_max-j+1) + (i_max-i)*j_max;
end

% function symk = sympmap(i, j, i_max, j_max)
%     symk = (j_max-j+1) + (i_max-i)*j_max;
% end