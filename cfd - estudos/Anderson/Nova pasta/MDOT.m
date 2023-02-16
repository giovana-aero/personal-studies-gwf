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
function continuity = MDOT(rho, u, rho_inf, u_inf, LVERT, DY)
% Check the validity of the numerical solution by confirming conservation of mass.
% The continuity equation in integral form is applied to the converged (steady) flow field and the rate
% of mass inflow across the entrance of the computational domain is compared to the rate of mass outflow across the exit boundary.
% The solution is considered valid when the deviation between the mass flow rate at the entrance and exit is less than 1 percent.
m_dot_in = - rho_inf*u_inf*(LVERT - DY) - DY*rho_inf*u_inf/2; % No need to integrate on the inflow boundary because
% both rho and u are constant except at the leading edge.
y = 0:DY:LVERT;
IMAX = size(rho,2);
m_dot_out = trapz(y, rho(:,IMAX).*u(:,IMAX));
deviation = abs((m_dot_in + m_dot_out)/m_dot_in)*100
if (deviation < 1)
    continuity = true;
else
    continuity = false;
end
