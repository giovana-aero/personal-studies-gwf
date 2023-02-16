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
function q_x = QX(T, k, DX, call_case)
% Calculate the x-component of the heat flux vector.
[JMAX, IMAX] = size(T);
dT_dx = zeros(JMAX, IMAX);
% Calculate the derivative of temperature wrt x:
% For interior points we have two cases. At the boundaries there is only one possibility: a forward difference
% when i = 1 and a rearward difference when i = IMAX
if (strcmp(call_case, 'Predict_E'))
    for i = 2:IMAX
        for j = 1:JMAX
            dT_dx(j,i) = (T(j,i) - T(j,i-1))/DX; % Rearward difference (opposite to predictor)
        end
    end
    dT_dx(:,1) = (T(:,2) - T(:,1))/DX; % Forward at i = 1
elseif (strcmp(call_case, 'Correct_E'))
    for i = 1:IMAX-1
        for j = 1:JMAX
            dT_dx(j,i) = (T(j,i+1) - T(j,i))/DX; % Forward difference (opposite to corrector)
        end
    end
    dT_dx(:,IMAX) = (T(:,IMAX) - T(:,IMAX-1))/DX; % Rearward at i = IMAX
else
    error('Undefined call case.')
end
% Then, following Fourier's law for heat conduction we have:
q_x = - k.*dT_dx;
