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
function [E1, E2, E3, E5] = Primitive2E(rho, u, p, v, T, mu, lambda, k, c_v, DX, DY, call_case)
% Calculate the flux vector E from the primitive flow field variables.
E1 = rho.*u; % From continuity
% To compute E2 first we need to calculate the x-direction of the normal stress (tau_xx)
tau_xx = TAUXX(u, v, lambda, mu, DX, DY, call_case)
% Then:
E2 = rho.*u.^2 + p - tau_xx; % From x-momentum
% For E3 we need the xy component of the shear stress (tau_xy)
tau_xy = TAUXY(u, v, mu, DX, DY, call_case)
E3 = rho.*u.*v - tau_xy; % From y-momentum
% Finally, for E5 we also need the x-component of the heat flux (q_x)
q_x = QX(T, k, DX, call_case)
E5 = (rho.*(c_v*T + (u.^2 + v.^2)/2) + p).*u - u.*tau_xx - v.*tau_xy + q_x; % From energy
