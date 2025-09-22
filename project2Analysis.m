%% Project 1 Analysis
% ------------------------------------------------------------------------------
% Adam Buchalter
%
% This is a script for:
% - Estimating and visualizing convergence of truncation error for the 2-D
%   cartesian slab multigroup neutron transport problem
% - Visualizing flux line profiles
% ------------------------------------------------------------------------------
clear all; close all; clc;

%% User Input
% ------------------------------------------------------------------------------

% Params
meshSizes = [129, 65, 33];
myCWD = pwd;

% Quick Maths
meshSizes = sort(meshSizes);
minN = meshSizes(1);
meshes = size(meshSizes, 2);
assert(meshes>=3, 'Need at Least 3 Grids')
midline = (minN+1)/2;

resFolder = 'results\\project2';

% ------------------------------------------------------------------------------

%% Grid Convergence Study
% ------------------------------------------------------------------------------

critdims = [];

for i = 1:meshes
    N = meshSizes(i);

    % File Info
    subfolder = 'results\\project2\\'+string(N)+'x'+string(N);
    resultmat = 'myTable.mat';
    resultOut = fullfile(myCWD,subfolder,resultmat);

    s = load(resultOut);
    dataTable = s.myTable;
    dataTable = sortrows(dataTable, 'keff');

    % num2str(s.myTable.keff,'%.5f') % Report to Requested Accuracy

    numRows = height(dataTable);
    fuelSizes = [];
    keffs = [];
    for j = 1:numRows
        fuelSize = dataTable.fuelSize(j);
        modThick = dataTable.modThick(j);
        % size = dataTable.size(j);
        keff = dataTable.keff(j);
        if modThick == 0
            fuelSizes = [fuelSizes, fuelSize];
            keffs = [keffs, keff];
        end
    end

    coefficients = polyfit(fuelSizes, keffs, 2);
    xFit = linspace(min(fuelSizes), max(fuelSizes));
    yFit = polyval(coefficients , xFit);
    objFun = coefficients;
    objFun(end) = objFun(end) - 1;
    r = min(roots(objFun));
    figure(i);
    % Plot flux surface
    plot(fuelSizes, keffs, '.', 'MarkerSize', 15);
    hold on;
    plot(xFit, yFit, 'k--'); % Plot Fit
    hold on;
    yl = yline(1,'--','$k_{eff}=1$','interpreter','latex');
    yl.LabelHorizontalAlignment = 'center';
    yl.Color = [.90 0 0];
    hold on;
    xl = xline(r,'--',[num2str(r,'%.4f'),' cm'],'interpreter','latex');
    xl.Color = [0 0 .90];
    hold off;
    ylabel('Predicted Neutron Multiplication Factor $k_{eff}$','interpreter','latex');
    xlabel('Fuel Slab Size [cm]','interpreter','latex');
    title(string(N)+'x'+string(N)+' Criticality Search','interpreter','latex');
    filename = 'critSearch.jpg';
    saveas(figure(i),fullfile(myCWD,subfolder,filename));

    critdims = [critdims, r];
end

% Init
f = zeros(meshes, 1);

% Evaluate Order of Convergence Directly
for i = 1:meshes
    f(i, 1) = critdims(i);
end
r = 2;
p2 = log(norm(f(1,1)-f(2,1))/norm(f(2,1)-f(3,1)))/log(r);
fprintf('Order of Grid Convergence by Critical Dimension: %g\n', p2);

% Demonstrate Order of Convergence Graphically
h = 1./meshSizes;
% Line of Fit
coefficients = polyfit(h, critdims, 1);
xFit = linspace(0, max(h));
yFit = polyval(coefficients , xFit);
extrapCritDim = coefficients(end);
% Inspect Closeness
figure(i+1);
plot(h, critdims, '.', 'MarkerSize', 15); % Plot Real Data
hold on;
plot(xFit, yFit, 'k--'); % Plot Fit
% hold on;
% yl = yline(1,'--','Analytical $k_{eff}$','interpreter','latex');
% yl.LabelHorizontalAlignment = 'center';
% yl.Color = [.90 0 0];
hold on;
yl = yline(extrapCritDim,'--',['Extrapolated: ',num2str(extrapCritDim,'%.2f'),' cm'],'interpreter','latex');
yl.LabelHorizontalAlignment = 'center';
yl.Color = [0 0 .90];
hold off;
ylabel('Calculated Critical Dimension [cm]','interpreter','latex');
xlabel('Relative Grid Spacing $h$','interpreter','latex');
title('First-Order Convergence in Critical Dimension','interpreter','latex');
set ( gca, 'XDir', 'reverse' )
filename = 'GridConvergence_critDim.jpg';
saveas(figure(i+1),fullfile(myCWD,resFolder,filename));

% ------------------------------------------------------------------------------

%% Reflector Savings Study
% ------------------------------------------------------------------------------

N = 129; % "Sufficient" discretization for converged solution (only accurate on order of 1 cm)

% File Info
subfolder = 'results\\project2\\'+string(N)+'x'+string(N);
resultmat = 'myTable.mat';
resultOut = fullfile(myCWD,subfolder,resultmat);

% plotsfolder = [num2str(fuelDim),'_',num2str(modDim),'_',num2str(size)];
% plotOut = fullfile(myCWD,subfolder,plotsfolder);

s = load(resultOut);
dataTable = s.myTable;
dataTable = sortrows(dataTable, 'fuelSize');

% num2str(s.myTable.keff,'%.5f') % Report to Requested Accuracy

numRows = height(dataTable);
fuelSizes = [];
modThicks = [];
keffs = [];
for i = 1:numRows
    fuelSize = dataTable.fuelSize(i);
    modThick = dataTable.modThick(i);
    % size = dataTable.size(i);
    keff = dataTable.keff(i);

    fuelSizes = [fuelSizes, fuelSize];
    modThicks = [modThicks, modThick];
    keffs = [keffs, keff];
end

myOrder = 3;

figure(10);
% Plot reflector optimization surface
scatter3(fuelSizes,modThicks,keffs);
hold on;
p = polyfitn([fuelSizes(:), modThicks(:)],keffs,myOrder);
myX = linspace(min(fuelSizes), max(fuelSizes));
myY = linspace(min(modThicks), max(modThicks));
[xg, yg] = meshgrid(myX,myY);
zg = polyvaln(p,[xg(:),yg(:)]);
surf(xg,yg,reshape(zg,size(xg)));
hold off;
ylabel('Moderator Thickness [cm]','interpreter','latex');
xlabel('Fuel Slab Dimension [cm]','interpreter','latex');
zlabel('Predicted Neutron Multiplication Factor $k_{eff}$','interpreter','latex')
title('Optimization of Reflector Savings');
filename = 'modelTraining.jpg';
saveas(figure(10),fullfile(myCWD,resFolder,filename));

% Find intersection with keff=1 plane
syms X1 X2
expression = polyn2sym(p);
eqn = expression == 1;
myModeratorVals = [];
myFuelVals = [];
for i = 1:length(myY)
    thisY = myY(i);
    eqn1 = subs(eqn, X2, thisY);
    sol_expr = solve(eqn1, X1);
    for j = 1:length(sol_expr)
        thisRoot = double(sol_expr(j));
        if isreal(thisRoot) && thisRoot <= max(fuelSizes) && thisRoot >= min(fuelSizes)
            myFuelVals = [myFuelVals, real(thisRoot)];
            myModeratorVals = [myModeratorVals, thisY];
        end
    end
end

figure(11);
plot(myModeratorVals, myFuelVals);
hold on;
[min_value, min_index] = min(myFuelVals);
minX = myModeratorVals(min_index);
xl = xline(minX,'--',['Optimal Thickness: ',num2str(minX,'%.2f'),' cm'],'interpreter','latex');
xl.Color = [0 0 .90];
hold off;
xlabel('Moderator Thickness [cm]','interpreter','latex');
ylabel('Fuel Slab Dimension [cm]','interpreter','latex');
title('Reflector Savings in a Critical System','interpreter','latex');
filename = 'reflectorSavings.jpg';
saveas(figure(11),fullfile(myCWD,resFolder,filename));

% ------------------------------------------------------------------------------