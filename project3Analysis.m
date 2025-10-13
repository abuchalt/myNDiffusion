%% Project 3 Analysis
% ------------------------------------------------------------------------------
% Adam Buchalter
%
% This is a script for:
% - Estimating, visualizing, and comparing convergence of truncation error
%   between explicit and implicit formulations of the 2-D cartesian slab
%   heat diffusion problem in steady-state and transient
% - Comparing computational efficiency (CPU-time) among the explicit and
%   implicit formulations
% - Visualizing temperature line profiles
% ------------------------------------------------------------------------------
clear all; close all; clc;

%% User Inputs
% ------------------------------------------------------------------------------

myCWD = pwd;
expFolder = fullfile(myCWD,'results\\project3exp');
impFolder = fullfile(myCWD,'results\\project3imp');
myOutFolder = fullfile(myCWD,'results\\project3');
mkdir(myOutFolder);

mySize = 12.5; % Domain Size [cm]

%% Explicit
% ------------------------------------------------------------------------------

% Spatial Convergence
% ------------------------------------------------------------------------------
meshSizes = [17, 33, 65];
timeStep = 0.005;

% Quick Maths
meshSizes = sort(meshSizes);
minN = meshSizes(1);
meshes = size(meshSizes, 2);
assert(meshes>=3, 'Need at Least 3 Grids')
midline = (minN+1)/2;

% Init
TPages = zeros(minN,minN,meshes);
PeakT = zeros(meshes, 1);

for i = 1:meshes
    N = meshSizes(i);

    % File Info
    folderName = string(N)+'x'+string(N)+'\\'+num2str(timeStep)+'dt_'+num2str(mySize)+'cm';
    Tmat = 'Tplot.mat';
    resultOut = fullfile(expFolder,folderName,Tmat);

    s = load(resultOut);
    Tplot = s.Tref2Plot;

    r = (N-1)/(minN-1);
    assert(floor(r)==r, 'Grid Refinement Ratio is Non-Integer');
    TPages(1:minN,1:minN,i) = Tplot(1:r:N,1:r:N);
    PeakT(i) = max(Tplot(:));
end

% Estimate Order of Convergence by peak Temp

% Evaluate Order of Convergence Directly
r = 2;
p2 = log(norm(PeakT(1,1)-PeakT(2,1))/norm(PeakT(2,1)-PeakT(3,1)))/log(r);
fprintf('Order of Grid Convergence by Peak Temperature: %g\n', p2);

% Demonstrate Order of Convergence Graphically
h = 1./meshSizes;
% Line of Fit
coefficients = polyfit(h, PeakT, 1);
xFit = linspace(0, max(h));
yFit = polyval(coefficients , xFit);
% Inspect Closeness
figure(1);
plot(h, PeakT, '.', 'MarkerSize', 15); % Plot Real Data
hold on;
plot(xFit, yFit, 'k--'); % Plot Fit
hold on;
extrapPeakT = polyval(coefficients,0);
yl = yline(extrapPeakT,'--',['Extrapolated: ',num2str(extrapPeakT,'%.2f'),' K'],'interpreter','latex');
yl.LabelHorizontalAlignment = 'center';
yl.LabelVerticalAlignment = 'bottom';
yl.Color = [.90 0 0];
hold off;
ylabel('Peak Slab Temperature [K]','interpreter','latex');
xlabel('Relative Grid Spacing $h$','interpreter','latex');
title('Estimation of Order of Spatial Convergence by Peak Steady-State Temperature','interpreter','latex');
set(gca, 'XDir', 'reverse');
filename = 'ExplicitSpatialConv.jpg';
saveas(figure(1),fullfile(myOutFolder,filename));

% Temporal Convergence
% ------------------------------------------------------------------------------
meshSize = 65;
timeSteps = [0.00125, 0.0025, 0.005];
T_infty = 18;
h = 0.1;

% Quick Maths
timeSteps = sort(timeSteps);
mindt = timeSteps(1);
meshes = size(timeSteps, 2);
maxdt = timeSteps(meshes);
assert(meshes>=3, 'Need at Least 3 Grids')

folderName = string(meshSize)+'x'+string(meshSize)+'\\'+num2str(mindt)+'dt_'+num2str(mySize)+'cm\\'+num2str(T_infty)+'T_'+num2str(h)+'h';
Tvec = 'transientT.mat';
resultOut = fullfile(expFolder,folderName,Tvec);

s = load(resultOut);
peakT = s.peakT;
minN = floor(size(peakT,2)/20);

% Init
TPages = zeros(minN,meshes);

for i = 1:meshes
    dt = timeSteps(i);

    % File Info
    folderName = string(meshSize)+'x'+string(meshSize)+'\\'+num2str(dt)+'dt_'+num2str(mySize)+'cm\\'+num2str(T_infty)+'T_'+num2str(h)+'h';
    Tvec = 'transientT.mat';
    TempOut = fullfile(expFolder,folderName,Tvec);
    % timevec = 'timeVec.mat';
    % TimeOut = fullfile(expFolder,folderName,timevec);

    s = load(TempOut);
    peakT = s.peakT;
    % s = load(TimeOut);
    % timeVec = s.timeVec;

    r = maxdt/dt;
    assert(floor(r)==r, 'Grid Refinement Ratio is Non-Integer');
    TPages(1:minN,i) = peakT(1:r:r*minN);
    % TPages(1:minN,1:minN,i) = Tplot(1:r:N,1:r:N);
    % PeakT(i) = max(Tplot(:));
end

figure(2);
for i = 1:meshes
    myErrContour = TPages(:,i)-TPages(:,1);
    errList(i) = sum(myErrContour,"all")/minN;
    plot(myErrContour);
    hold on;
end
hold off;
ylabel('Relative Error','interpreter','latex');
xlabel('Relative Time','interpreter','latex');
title('Global Error Reduction with $\Delta t$ Refinement','interpreter','latex');
legend('0.00125 s', '0.0025 s', '0.005 s');
filename = 'ExpGlobalErrorReduction.jpg';
saveas(figure(2),fullfile(myOutFolder,filename));

% Evaluate Order of Convergence Directly
r = 2;
p2 = log(norm(errList(3)-errList(2))/norm(errList(2)-errList(1)))/log(r);
fprintf('Order of Grid Convergence by Transient Global Error: %g\n', p2);

% Demonstrate Order of Convergence Graphically
h = timeSteps;
% Line of Fit
coefficients = polyfit(h, errList, 1);
xFit = linspace(0, max(h));
yFit = polyval(coefficients , xFit);
% Inspect Closeness
figure(3);
plot(h, errList*1E4, '.', 'MarkerSize', 15); % Plot Real Data
hold on;
plot(xFit, yFit*1E4, 'k--'); % Plot Fit
% extrapPeakT = polyval(coefficients,0);
% yl = yline(extrapPeakT,'--',['Extrapolated: ',num2str(extrapPeakT,'%.2f'),' K'],'interpreter','latex');
% yl.LabelHorizontalAlignment = 'center';
% yl.LabelVerticalAlignment = 'bottom';
% yl.Color = [.90 0 0];
hold off;
ylabel('Relative Global Error $\times 10^{-4}$','interpreter','latex');
xlabel('Relative Time Spacing $\Delta t$','interpreter','latex');
title('Estimate Temporal Convergence Order by Global Error Reduction in Transient','interpreter','latex');
set(gca, 'XDir', 'reverse');
filename = 'ExplicitTemporalConv.jpg';
saveas(figure(3),fullfile(myOutFolder,filename));


%% Implicit
% ------------------------------------------------------------------------------

% Spatial Convergence
% ------------------------------------------------------------------------------
meshSizes = [17, 33, 65];
timeStep = 0.1;

% Quick Maths
meshSizes = sort(meshSizes);
minN = meshSizes(1);
meshes = size(meshSizes, 2);
assert(meshes>=3, 'Need at Least 3 Grids')
midline = (minN+1)/2;

% Init
TPages = zeros(minN,minN,meshes);
PeakT = zeros(meshes, 1);

for i = 1:meshes
    N = meshSizes(i);

    % File Info
    folderName = string(N)+'x'+string(N)+'\\'+num2str(timeStep)+'dt_'+num2str(mySize)+'cm';
    Tmat = 'Tplot.mat';
    resultOut = fullfile(impFolder,folderName,Tmat);

    s = load(resultOut);
    Tplot = s.Tref2Plot;

    r = (N-1)/(minN-1);
    assert(floor(r)==r, 'Grid Refinement Ratio is Non-Integer');
    TPages(1:minN,1:minN,i) = Tplot(1:r:N,1:r:N);
    PeakT(i) = max(Tplot(:));
end

% Estimate Order of Convergence by peak Temp

% Evaluate Order of Convergence Directly
r = 2;
p2 = log(norm(PeakT(1,1)-PeakT(2,1))/norm(PeakT(2,1)-PeakT(3,1)))/log(r);
fprintf('Order of Grid Convergence by Peak Temperature: %g\n', p2);

% Demonstrate Order of Convergence Graphically
h = 1./meshSizes;
% Line of Fit
coefficients = polyfit(h, PeakT, 1);
xFit = linspace(0, max(h));
yFit = polyval(coefficients , xFit);
% Inspect Closeness
figure(4);
plot(h, PeakT, '.', 'MarkerSize', 15); % Plot Real Data
hold on;
plot(xFit, yFit, 'k--'); % Plot Fit
hold on;
extrapPeakT = polyval(coefficients,0);
yl = yline(extrapPeakT,'--',['Extrapolated: ',num2str(extrapPeakT,'%.2f'),' K'],'interpreter','latex');
yl.LabelHorizontalAlignment = 'center';
yl.LabelVerticalAlignment = 'bottom';
yl.Color = [.90 0 0];
hold off;
ylabel('Peak Slab Temperature [K]','interpreter','latex');
xlabel('Relative Grid Spacing $h$','interpreter','latex');
title('Estimation of Order of Spatial Convergence by Peak Steady-State Temperature','interpreter','latex');
set(gca, 'XDir', 'reverse');
filename = 'ImplicitSpatialConv.jpg';
saveas(figure(4),fullfile(myOutFolder,filename));

% Temporal Convergence
% ------------------------------------------------------------------------------
meshSize = 65;
timeSteps = [1, 2, 4];
T_infty = 18;
h = 0.1;

% Quick Maths
timeSteps = sort(timeSteps);
mindt = timeSteps(1);
meshes = size(timeSteps, 2);
maxdt = timeSteps(meshes);
assert(meshes>=3, 'Need at Least 3 Grids')

folderName = string(meshSize)+'x'+string(meshSize)+'\\'+num2str(mindt)+'dt_'+num2str(mySize)+'cm\\'+num2str(T_infty)+'T_'+num2str(h)+'h';
Tvec = 'transientT.mat';
resultOut = fullfile(impFolder,folderName,Tvec);

s = load(resultOut);
peakT = s.peakT;
minN = floor(size(peakT,2)/20);

% Init
TPages = zeros(minN,meshes);

for i = 1:meshes
    dt = timeSteps(i);

    % File Info
    folderName = string(meshSize)+'x'+string(meshSize)+'\\'+num2str(dt)+'dt_'+num2str(mySize)+'cm\\'+num2str(T_infty)+'T_'+num2str(h)+'h';
    Tvec = 'transientT.mat';
    TempOut = fullfile(impFolder,folderName,Tvec);
    % timevec = 'timeVec.mat';
    % TimeOut = fullfile(impFolder,folderName,timevec);

    s = load(TempOut);
    peakT = s.peakT;
    % s = load(TimeOut);
    % timeVec = s.timeVec;

    r = maxdt/dt;
    assert(floor(r)==r, 'Grid Refinement Ratio is Non-Integer');
    TPages(1:minN,i) = peakT(1:r:r*minN);
    % TPages(1:minN,1:minN,i) = Tplot(1:r:N,1:r:N);
    % PeakT(i) = max(Tplot(:));
end

figure(5);
for i = 1:meshes
    myErrContour = TPages(:,i)-TPages(:,1);
    errList(i) = sum(myErrContour,"all")/minN;
    plot(myErrContour);
    hold on;
end
hold off;
ylabel('Relative Error','interpreter','latex');
xlabel('Relative Time','interpreter','latex');
title('Global Error Reduction with $\Delta t$ Refinement','interpreter','latex');
legend('0.00125 s', '0.0025 s', '0.005 s');
filename = 'ImpGlobalErrorReduction.jpg';
saveas(figure(5),fullfile(myOutFolder,filename));

% Evaluate Order of Convergence Directly
r = 2;
p2 = log(norm(errList(3)-errList(2))/norm(errList(2)-errList(1)))/log(r);
fprintf('Order of Grid Convergence by Transient Global Error: %g\n', p2);

% Demonstrate Order of Convergence Graphically
h = timeSteps;
% Line of Fit
coefficients = polyfit(h, errList, 1);
xFit = linspace(0, max(h));
yFit = polyval(coefficients , xFit);
% Inspect Closeness
figure(6);
plot(h, errList*1E3, '.', 'MarkerSize', 15); % Plot Real Data
hold on;
plot(xFit, yFit*1E3, 'k--'); % Plot Fit
% extrapPeakT = polyval(coefficients,0);
% yl = yline(extrapPeakT,'--',['Extrapolated: ',num2str(extrapPeakT,'%.2f'),' K'],'interpreter','latex');
% yl.LabelHorizontalAlignment = 'center';
% yl.LabelVerticalAlignment = 'bottom';
% yl.Color = [.90 0 0];
hold off;
ylabel('Relative Global Error $\times 10^{-3}$','interpreter','latex');
xlabel('Relative Time Spacing $\Delta t$','interpreter','latex');
title('Estimate Temporal Convergence Order by Global Error Reduction in Transient','interpreter','latex');
set(gca, 'XDir', 'reverse');
filename = 'ImplicitTemporalConv.jpg';
saveas(figure(6),fullfile(myOutFolder,filename));