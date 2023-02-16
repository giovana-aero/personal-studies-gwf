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
function tau_xy = TAUXY(u, v, mu, DX, DY, call_case)
% Calculate the xy (or yx) component of the shear stress.
[JMAX, IMAX] = size(u);
du_dy = zeros(JMAX, IMAX);
dv_dx = zeros(JMAX, IMAX);
% Calculate the derivative of u wrt y and the derivative of v wrt x:
% For interior points we have a total of four cases combining the different finite difference schemes used to calculate du_dy and dv_dx.
% A the boundaries there is only one possibility: a forward difference when i,j = 1 and a rearward difference when i,j = IMAX,JMAX.
if (strcmp(call_case, 'Predict_E') || strcmp(call_case, 'Correct_E'))
    % The calculation of the derivative of u wrt y is the same for these two cases:
    for i = 1:IMAX
        for j = 2:JMAX-1
            du_dy(j,i) = (u(j+1,i) - u(j-1,i))/(2*DY); % Central difference
        end
    end
    du_dy(1,:) = (u(2,:) - u(1,:))/DY; % Forward at j = 1
    du_dy(JMAX,:) = (u(JMAX,:) - u(JMAX-1,:))/DY; % Rearward at j = JMAX

    % The calculation of the derivative of v wrt x is different for these two cases:
    if (strcmp(call_case, 'Predict_E'))
        % Rearward difference except at i = 1
        for i = 2:IMAX
            for j = 1:JMAX
                dv_dx(j,i) = (v(j,i) - v(j,i-1))/DX;
            end
        end
        dv_dx(:,1) = (v(:,2) - v(:,1))/DX; % Forward at i = 1
    else
        % Forward difference except at i = IMAX
        for i = 1:IMAX-1
            for j = 1:JMAX
                dv_dx(j,i) = (v(j,i+1) - v(j,i))/DX;
            end
        end
        dv_dx(:,IMAX) = (v(:,IMAX) - v(:,IMAX-1))/DX; % Rearward at i = IMAX
    end
elseif (strcmp(call_case, 'Predict_F') || strcmp(call_case, 'Correct_F'))
    % The calculation of the derivative of v wrt x is the same for these two cases:
    for i = 2:IMAX-1
        for j = 1:JMAX
            dv_dx(j,i) = (v(j,i+1) - v(j,i-1))/(2*DX); % Central difference
        end
    end
    dv_dx(:,1) = (v(:,2) - v(:,1))/DX; % Forward at i = 1
    dv_dx(:,IMAX) = (v(:,IMAX) - v(:,IMAX-1))/DX; % Rearward at i = IMAX

    % The calculation of the derivative of u wrt y is different for these two cases:
    if (strcmp(call_case, 'Predict_F'))
        % Rearward difference except at j = 1
        for i = 1:IMAX
            for j = 2:JMAX
                du_dy(j,i) = (u(j,i) - u(j-1,i))/DY;
            end
        end
        du_dy(1,:) = (u(2,:) - u(1,:))/DY; % Forward at j = 1
    else
        % Forward difference except at j = JMAX
        for i = 1:IMAX
            for j = 1:JMAX-1
                du_dy(j,i) = (u(j+1,i) - u(j,i))/DY;
            end
        end
        du_dy(JMAX,:) = (u(JMAX,:) - u(JMAX-1,:))/DY; % Rearward at j = JMAX
    end
else
    error('Undefined call case.')
end
% Finally, from the shear stress definition (assuming a Newtonian fluid) we have:
tau_xy = mu.*(du_dy + dv_dx);
