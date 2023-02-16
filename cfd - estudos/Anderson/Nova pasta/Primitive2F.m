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
function [F1, F2, F3, F5] = Primitive2F(rho, u, p, v, T, mu, lambda, k, c_v, DX, DY, call_case)
% Calculate the flux vector F from the primitive flow field variables.
F1 = rho.*v; % From continuity
% For F2 we need the yx component of the shear stress (tau_yx)
tau_yx = TAUXY(u, v, mu, DX, DY, call_case) % Atteding to the fact that tau_yx = tau_xy, the same function can be used to compute both shear stresses.
F2 = rho.*u.*v - tau_yx; % From x-momentum
% To compute F3 first we need to calculate the y-direction of the normal stress (tau_yy)
tau_yy = TAUYY(u, v, lambda, mu, DX, DY, call_case)
% Then:
F3 = rho.*v.^2 + p - tau_yy; % From y-momentum
% Finally, for F5 we also need the y-component of the heat flux (q_y)
q_y = QY(T, k, DY, call_case),disp(q_y(5,2))
F5 = (rho.*(c_v*T + (u.^2 + v.^2)/2) + p).*v - u.*tau_yx - v.*tau_yy + q_y;
