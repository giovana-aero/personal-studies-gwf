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
function [rho, u, v, p, T] = BC(rho, u, v, p, T, rho_inf, u_inf, p_inf, T_inf, T_w_T_inf, R, adiabatic_wall_flag)
% Apply boundary conditions (constant-temperature wall case).
[JMAX, IMAX] = size(rho);
% We have four different cases:
% Case 1: Leading edge --> No-slip condition is enforced on velocity and the rest of the primitive flow field
% variables are assumed to be equal to their respective freestream values:
% u(1,1) = 0 and v(1,1) = 0 (implicit in the definition of u and v)
T(1,1) = T_inf;
p(1,1) = p_inf;
rho(1,1) = rho_inf;
% Case 2: Inflow / Upper boundary (not leading edge) --> y component of velocity (v) is assumed equal to zero
% and the rest of flow field properties are assumed to be their respective freestream values:
% For the inflow boundary:
u(2:JMAX,1) = u_inf;
% v(2:JMAX,1) = 0 is implicit in the definition of v
p(2:JMAX,1) = p_inf;
T(2:JMAX,1) = T_inf;
rho(2:JMAX,1) = rho_inf;
% For the upper boundary:
u(JMAX,2:IMAX) = u_inf;
% v(JMAX,2:IMAX) = 0 is implicit in the definition of v
p(JMAX,2:IMAX) = p_inf;
T(JMAX,2:IMAX) = T_inf;
rho(JMAX,2:IMAX) = rho_inf;
% Case 4: Outflow (not surface or JMAX) --> Velocity, pressure and
% temperature are calculated based on an extrapolation from the two adjacent interior points.
% Density is computed from the equation of state:
u(2:JMAX-1,IMAX) = 2*u(2:JMAX-1,IMAX-1) - u(2:JMAX-1,IMAX-2);
v(2:JMAX-1,IMAX) = 2*v(2:JMAX-1,IMAX-1) - v(2:JMAX-1,IMAX-2);
p(2:JMAX-1,IMAX) = 2*p(2:JMAX-1,IMAX-1) - p(2:JMAX-1,IMAX-2);
T(2:JMAX-1,IMAX) = 2*T(2:JMAX-1,IMAX-1) - T(2:JMAX-1,IMAX-2);
rho(2:JMAX-1,IMAX) = p(2:JMAX-1,IMAX)./(R*T(2:JMAX-1,IMAX));
% Note: Case 4 must be executed before case 3 because the pressure
% extrapolation at IMAX needs the values of the outflow boundary
% Case 3: Surface (not leading edge) --> No-slip condition is specified on velocity.
% Temperature is either assumed to be equal to the wall temperature value
% or the adiabatic wall boundary condition is applied.
% Pressure is calculated by extrapolating from the values at the two points adjacent to the surface.
% Density is computed from the equation of state:
% u(1,2:IMAX) = 0 and v(1,2:IMAX) = 0 (implicit in the definition of u and v)
if (adiabatic_wall_flag) % Then an adiabatic wall is assumed
    T(1,2:IMAX) = T(2,2:IMAX); % This condition results from enforcing that the normal temperature gradient pointing into to the wall is zero.
    % In the present case this condition is traduced in dT_dy = 0 at the wall.
else
    T(1,2:IMAX) = T_w_T_inf*T_inf;
end
p(1,2:IMAX) = 2*p(2,2:IMAX) - p(3,2:IMAX);
rho(1,2:IMAX) = p(1,2:IMAX)./(R*T(1,2:IMAX));
