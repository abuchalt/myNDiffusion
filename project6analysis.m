%% 2D Coupled SMR Results Visualizer
% ------------------------------------------------------------------------------
% Adam Buchalter
%
% desc
% ------------------------------------------------------------------------------
clear all; close all; clc;
myCWD = pwd;

%% User Input
% ------------------------------------------------------------------------------

% File Info
% subfolder='results\\project6\\shellFuel_TM557';
% subfolder='results\\project6\\checkeredFuel_TM557';
% subfolder='results\\project6\\uniformFuel_TM557';
mkdir(fullfile(myCWD,subfolder));

%% Params
% ------------------------------------------------------------------------------

% Physical params
totPwr = 250*1E6; % Thermal Output [MW_th -> W_th]
% NbyN = 17; % 17x17 Fuel Assembly
fuelLength = 78.74*2.54; % Active Height [in -> cm]
pitch_Assy = 8.466*2.54; % Assembly Pitch [in -> cm]
totLinPwr = totPwr/fuelLength; % Total Linear Heat Generation [W/cm]

% Solved params
S = load(fullfile(subfolder,'myVals.mat'));
TM_nom = S.TM_nom;
keff_iter = S.keff_iter;
myBOR = S.myBOR;

phi = load(fullfile(subfolder,'phi.mat')).phi;
TPlot = load(fullfile(subfolder,'TPlot.mat')).Tref2Plot;
w = load(fullfile(subfolder,'w.mat')).w;
xPlot = load(fullfile(subfolder,'xPlot.mat')).fullx;
yPlot = load(fullfile(subfolder,'yPlot.mat')).fully;

% Extracted Params
i_max = (size(xPlot, 1)+1)/2;
j_max = (size(xPlot, 2)+1)/2;
Deltax = xPlot(2) - xPlot(1);
Deltay = Deltax;

q3prime = phi(1:i_max*j_max)+phi(1+i_max*j_max:2*i_max*j_max);
q3prime = (totLinPwr/4)*q3prime/(sum(q3prime,'all')*Deltax*Deltay);

% Thermal Flux
phiPlot = reshape(phi(1+i_max*j_max:2*i_max*j_max), i_max, j_max);
phiref1 = rot90(phiPlot,3);
phiref1 = phiref1(:,1:end-1);
phiref1Plot = horzcat(phiref1, phiPlot);
phiref2 = rot90(phiref1Plot, 2);
phiref2 = phiref2(1:end-1,:);
phiref2Plot = vertcat(phiref2, phiref1Plot);
% Fast Flux
FphiPlot = reshape(phi(1:i_max*j_max), i_max, j_max);
Fphiref1 = rot90(FphiPlot,3);
Fphiref1 = Fphiref1(:,1:end-1);
Fphiref1Plot = horzcat(Fphiref1, FphiPlot);
Fphiref2 = rot90(Fphiref1Plot, 2);
Fphiref2 = Fphiref2(1:end-1,:);
Fphiref2Plot = vertcat(Fphiref2, Fphiref1Plot);

% Volumetric Heat Rate
q3primePlot = reshape(q3prime, i_max, j_max);
q3primeref1 = rot90(q3primePlot,3);
q3primeref1 = q3primeref1(:,1:end-1);
q3primeref1Plot = horzcat(q3primeref1, q3primePlot);
q3primeref2 = rot90(q3primeref1Plot, 2);
q3primeref2 = q3primeref2(1:end-1,:);
q3primeref2Plot = vertcat(q3primeref2, q3primeref1Plot);

% Temperature
TPlot;

% Coolant Velocity
wPlot = reshape(w, i_max, j_max);
wref1 = rot90(wPlot,3);
wref1 = wref1(:,1:end-1);
wref1Plot = horzcat(wref1, wPlot);
wref2 = rot90(wref1Plot, 2);
wref2 = wref2(1:end-1,:);
wref2Plot = vertcat(wref2, wref1Plot);

%% Plotting
% ------------------------------------------------------------------------------

figure(1);
% Plot centerline temp profile
plot(xPlot, TPlot(i_max,:))
% plot(xPlot, diag(TPlot))
ylabel('Temperature [K]');
xlabel('x [cm]');
title('Steady-State Centerline Temperature Profile');
dim = [.4 .3 .3 .3];
str = compose('TMO = %.0f',TM_nom);
annotation('textbox',dim,'String',str,'FitBoxToText','on');

figure(2);
hold on;
% Plot centerline flux profile
plot(xPlot, phiref2Plot(i_max,:), 'r');
hold on;
plot(xPlot, Fphiref2Plot(i_max,:), 'b');
ylabel('Normalized Flux');
xlabel('x [cm]');
title('Steady-State Centerline Flux Profiles');
dim = [.4 .3 .3 .3];
dim2 = [.4 .2 .3 .3];
str = compose('k_{eff} = %.6f',keff_iter);
str2 = compose('BOR = %.3f',myBOR);
annotation('textbox',dim,'String',str,'FitBoxToText','on');
annotation('textbox',dim2,'String',str2,'FitBoxToText','on');

figure(3);
% Plot centerline power profile
plot(xPlot, q3primeref2Plot(i_max,:))
ylabel('Volumetric Heat Rate [W/cc]');
xlabel('x [cm]');
title('Steady-State Centerline Power Profile');

figure(4);
% Plot centerline coolant velocity profile
plot(xPlot, wref2Plot(i_max,:))
ylabel('Coolant Velocity [cm/s]');
xlabel('x [cm]');
title('Steady-State Centerline Coolant Velocity Profile');

saveas(figure(1),fullfile(myCWD,subfolder,'tempProfile.jpg'));
saveas(figure(2),fullfile(myCWD,subfolder,'fluxProfile.jpg'));
saveas(figure(3),fullfile(myCWD,subfolder,'powerProfile.jpg'));
saveas(figure(4),fullfile(myCWD,subfolder,'velocityProfile.jpg'));