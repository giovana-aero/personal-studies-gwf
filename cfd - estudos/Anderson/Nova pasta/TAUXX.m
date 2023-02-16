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
function tau_xx = TAUXX(u, v, lambda, mu, DX, DY, call_case)
% Calculate the x-direction of the normal stress.
[JMAX, IMAX] = size(u);
du_dx = zeros(JMAX, IMAX);
dv_dy = zeros(JMAX, IMAX);
% Calculate the derivative of u wrt x and the derivative of v wrt y (divergence of 2D velocity):
% For interior points we have two cases for the derivative of u wrt x and a unique case for the derivative of v wrt y.
% At the boundaries there is only one possibility: a forward difference when i,j = 1 and a rearward difference when i,j = IMAX,JMAX.
% To compute du_dx:
if (strcmp(call_case, 'Predict_E'))
    for i = 2:IMAX
        for j = 1:JMAX
            du_dx(j,i) = (u(j,i) - u(j,i-1))/DX; % Rearward difference (opposite to predictor)
        end
    end
    du_dx(:,1) = (u(:,2) - u(:,1))/DX; % Forward at i = 1
elseif (strcmp(call_case, 'Correct_E'))
    for i = 1:IMAX-1
        for j = 1:JMAX
            du_dx(j,i) = (u(j,i+1) - u(j,i))/DX; % Forward difference (opposite to corrector)
        end
    end
    du_dx(:,IMAX) = (u(:,IMAX) - u(:,IMAX-1))/DX; % Rearward at i = IMAX
else
    error('Undefined call case.')
end
% To compute dv_dy:
for i = 1:IMAX
    for j = 2:JMAX-1
        dv_dy(j,i) = (v(j+1,i) - v(j-1,i))/(2*DY); % Central difference
    end
end
dv_dy(1,:) = (v(2,:) - v(1,:))/DY; % Forward at j = 1
dv_dy(JMAX,:) = (v(JMAX,:) - v(JMAX-1,:))/DY; % Rearward at j = JMAX
% Finally, from the normal stress definition (assuming a Newtonian fluid and the Stokes' hypothesis) we have:
tau_xx = lambda.*(du_dx + dv_dy) + 2*mu.*du_dx;
