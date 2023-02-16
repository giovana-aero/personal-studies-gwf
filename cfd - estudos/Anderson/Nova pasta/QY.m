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
function q_y = QY(T, k, DY, call_case)
% Calculate the y-component of the heat flux vector.
[JMAX, IMAX] = size(T);
dT_dy = zeros(JMAX, IMAX);
% Calculate the derivative of temperature wrt y:
% For interior points we have two cases. At the boundaries there is only one possibility: a forward difference
% when j = 1 and a rearward difference when j = JMAX
if (strcmp(call_case, 'Predict_F'))
    for i = 1:IMAX
        for j = 2:JMAX
            dT_dy(j,i) = (T(j,i) - T(j-1,i))/DY; % Rearward difference (opposite to predictor)
        end
    end
    dT_dy(1,:) = (T(2,:) - T(1,:))/DY; % Forward at j = 1
elseif (strcmp(call_case, 'Correct_F'))
    for i = 1:IMAX
        for j = 1:JMAX-1
            dT_dy(j,i) = (T(j+1,i) - T(j,i))/DY; % Forward difference (opposite to corrector)
        end
    end
    dT_dy(JMAX,:) = (T(JMAX,:) - T(JMAX-1,:))/DY; % Rearward at j = JMAX
else
    error('Undefined call case.')
end
% Then, following Fourier's law for heat conduction we have:
q_y = - k.*dT_dy;
