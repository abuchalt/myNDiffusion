Tplot = reshape(T, i_max, j_max);
Tref1 = rot90(Tplot,3);
Tref1 = Tref1(:,1:end-1);
Tref1Plot = horzcat(Tref1, Tplot);
Tref2 = rot90(Tref1Plot, 2);
Tref2 = Tref2(1:end-1,:);
Tref2Plot = vertcat(Tref2, Tref1Plot);
Tref2Plot = T_r*Tref2Plot; % Renormalize to Reference Temperature

figure(1);
% Plot temp surface
surf(fullx,fully,Tref2Plot);
ylabel('y [cm]');
xlabel('x [cm]');
zlabel('Temperature [K]')
title('Steady-State Solution for a Uniformly Heated Slab');