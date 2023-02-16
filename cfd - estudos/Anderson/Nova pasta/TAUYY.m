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
function tau_yy = TAUYY(u, v, lambda, mu, DX, DY, call_case)
% Calculate the y-direction of the normal stress.
[JMAX, IMAX] = size(v);
du_dx = zeros(JMAX, IMAX);
dv_dy = zeros(JMAX, IMAX);
% Calculate the derivative of u wrt x and the derivative of v wrt y (divergence of 2D velocity):
% For interior points we have two cases for the derivative of v wrt y and a unique case for the derivative of u wrt x.
% At the boundaries there is only one possibility: a forward difference when i,j = 1 and a rearward difference when i,j = IMAX,JMAX.
% To compute dv_dy:
if (strcmp(call_case, 'Predict_F'))
    for i = 1:IMAX
        for j = 2:JMAX
            dv_dy(j,i) = (v(j,i) - v(j-1,i))/DY; % Rearward difference (opposite to predictor)
        end
    end
    dv_dy(1,:) = (v(2,:) - v(1,:))/DY; % Forward at j = 1
elseif (strcmp(call_case, 'Correct_F'))
    for i = 1:IMAX
        for j = 1:JMAX-1
            dv_dy(j,i) = (v(j+1,i) - v(j,i))/DY; % Forward difference (opposite to corrector)
        end
    end
    dv_dy(JMAX,:) = (v(JMAX,:) - v(JMAX-1,:))/DY; % Rearward at j = JMAX
else
    error('Undefined call case.')
end
% To compute du_dx:
for i = 2:IMAX-1
    for j = 1:JMAX
        du_dx(j,i) = (u(j,i+1) - u(j,i-1))/(2*DX); % Central difference
    end
end
du_dx(:,1) = (u(:,2) - u(:,1))/DX; % Forward at i = 1
du_dx(:,IMAX) = (u(:,IMAX) - u(:,IMAX-1))/DX; % Rearward at i = IMAX
% Finally, from the normal stress definition (assuming a Newtonian fluid and the Stokes' hypothesis) we have:
tau_yy = lambda.*(du_dx + dv_dy) + 2*mu.*dv_dy;
