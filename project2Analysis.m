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
meshSizes = [257, 129, 65];
myCWD = pwd;

% Quick Maths
meshSizes = sort(meshSizes);
minN = meshSizes(1);
meshes = size(meshSizes, 2);
assert(meshes>=3, 'Need at Least 3 Grids')
midline = (minN+1)/2;

for i = 1:meshes
    N = meshSizes(i);

    % File Info
    subfolder = 'results\\project2\\'+string(N)+'x'+string(N);
    resultmat = 'myTable.mat';
    resultOut = fullfile(myCWD,subfolder,resultmat);

    % plotsfolder = [num2str(fuelDim),'_',num2str(modDim),'_',num2str(size)];
    % plotOut = fullfile(myCWD,subfolder,plotsfolder);

    s = load(resultOut);
    dataTable = s.myTable;
    dataTable = sortrows(dataTable, 'keff');

    num2str(s.myTable.keff,'%.5f') % Report to Requested Accuracy

    numRows = height(dataTable);
    fuelSizes = [];
    keffs = [];
    for i = 1:numRows
        fuelSize = dataTable.fuelSize(i);
        modThick = dataTable.modThick(i);
        size = dataTable.size(i);
        keff = dataTable.keff(i);
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
end