%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NUMERICAL SOLUTION OF SUPERSONIC FLOW OVER A FLAT PLATE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                         %
% Authors:                                                %
% Alberto Garcia (albertoomniaf1@gmail.com)               %
% Ivan Padilla (ivan.padilla6@gmail.com)                  %
% ------------------------------------------------------- %
% December 2013                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate plots for comparison with the figures in Anderson Chapter 10.
close all
clear all
clc
main_constant_T_wall
figure
plot(1:IMAX,p(1,:)/p_inf)
grid on
xlabel('Grid location along surface of plate')
ylabel('P/P_{\infty}')
title('Normalized surface pressure distributions')
y_norm_TE = y/x(IMAX)*sqrt(Re_L);
figure
plot(p(:,IMAX)/p_inf,y_norm_TE)
grid on
xlabel('P/P_{\infty} (At the trailing edge)')
ylabel('Normalized y distance')
title('Normalized pressure profiles')
figure
plot(T(:,IMAX)/T_inf,y_norm_TE)
grid on
xlabel('T/T_{\infty} (At the trailing edge)')
ylabel('Normalized y distance')
title('Normalized temperature profiles')
figure
plot(u(:,IMAX)/(M_inf*a_inf),y_norm_TE)
grid on
xlabel('u/u_{\infty} (At the trailing edge)')
ylabel('Normalized y distance')
title('Normalized velocity profiles')
figure
plot(M(:,IMAX),y_norm_TE)
grid on
xlabel('Local Mach number (At the trailing edge)')
ylabel('Normalized y distance')
title('Local Mach number profiles')
clear all
clc
main_adiabatic_wall
y_norm_TE = y/x(IMAX)*sqrt(Re_L);
figure(2)
hold on
plot(1:IMAX,p(1,:)/p_inf,'r')
legend('Constant T', 'Adiabatic')
figure(3)
hold on
plot(p(:,IMAX)/p_inf,y_norm_TE,'r')
legend('Constant T', 'Adiabatic')
figure(4)
hold on
plot(T(:,IMAX)/T_inf,y_norm_TE,'r')
legend('Constant T', 'Adiabatic')
figure(5)
hold on
plot(u(:,IMAX)/(M_inf*a_inf),y_norm_TE,'r')
legend('Constant T', 'Adiabatic')
figure(6)
hold on
plot(M(:,IMAX),y_norm_TE,'r')
legend('Constant T', 'Adiabatic')
