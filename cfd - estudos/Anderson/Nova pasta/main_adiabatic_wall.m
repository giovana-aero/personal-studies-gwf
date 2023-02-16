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
%%%%%%%%%%%%%
% Grid size %
%%%%%%%%%%%%%
IMAX = 70; % Number of grid points along the x-coordinate
JMAX = 70; % Number of grid points along the y-coordinate
MAXIT = 1e4; % Maximum number of time steps allowed before stopping the program
t = 1; % Time step counter
time = 0; % Physical time [s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Freestream conditions / Plate length %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ISA sea level conditions
M_inf = 4; % Freestream Mach number
a_inf = 340.28; % Freestream speed of sound [m/s]
p_inf = 101325; % Freestream pressure [Pa]
T_inf = 288.16; % Freestream temperature [K]
LHORI = 1e-5; % Flat plate length [m]
%%%%%%%%%%%%%%%%%%%%%%%%
% Constants definition %
%%%%%%%%%%%%%%%%%%%%%%%%
gamma = 1.4; % Ratio of specific heats for calorically perfect air
mu_0 = 1.7894e-5; % Reference (sea level) value for dynamic viscosity [kg/(m*s)]
T_0 = 288.16; % Reference (sea level) value for temperature [K]
Pr = 0.71; % Prandtl number for calorically perfect air
R = 287; % Specific gas constant for calorically perfect air [J/(kg*K)]
% Compute mu_inf by means of Sutherland's law by using the DYNVIS function
% mu_inf = DYNVIS(T_inf, mu_0, T_0); % Freestream dynamic viscosity [kg/(m*s)]
% Continue with the definition of other necessary variables
c_v = R/(gamma - 1); % Specific heat at constant volume for calorically perfect air [J/(kg*K)]
c_p = gamma*c_v; % Specific heat at constant pressure for calorically perfect air [J/(kg*K)]
rho_inf = p_inf/(R*T_inf); % Freestream density [kg/m^3]
Re_L = rho_inf*M_inf*a_inf*LHORI/DYNVIS(T_inf, mu_0, T_0); % Freestream Reynolds number according to the flat plate length
% e_inf = c_v*T_inf; % Freestream specific internal energy [J/kg]
% Compute k_inf by means of the constant Prandtl number assumption
% k_inf = THERMC(Pr, c_p, mu_inf); % Freestream thermal conductivity [W/(m*K)]
% Continue with the definition of other necessary variables
delta = 5*LHORI/sqrt(Re_L); % Initial assumption of the boundary layer thickness as predicted by a
% Blasius calculation at the trailind edge [m]
LVERT = 5*delta; % Assumed height for the computational domain [m]
% Calculate step sizes in x and y
DX = LHORI/(IMAX - 1); % Spacing between grid points along the x-coordinate [m]
DY = LVERT/(JMAX - 1); % Spacing between grid points along the y-coordinate [m]
x = 0:DX:LHORI; % x-coordinate of the domain
y = 0:DY:LVERT; % y-coordinate of the domain
CFL = 0.5; % Courant-Friedrichs-Lewy number (should be between 0.5 and 0.8)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flow field variables initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Only those variables that are needed in the calculation of the vectors U, E and F are initialized.
% The remaining flow properties can be obtained after as a postprocess to the solution.
p = ones(JMAX,IMAX)*p_inf;
rho = ones(JMAX,IMAX)*rho_inf;
T = ones(JMAX,IMAX)*T_inf;
u = ones(JMAX,IMAX)*M_inf*a_inf;
u(1,:) = 0; % No-slip boundary condition enforced at the surface
v = zeros(JMAX,IMAX);
mu = DYNVIS(T, mu_0, T_0);
lambda = - 2/3*mu; % Second viscosity (based on Stokes' hypothesis) [kg/(m*s)]
k = THERMC(Pr, c_p, mu);
% Other variables needed for intermediate calculations are defined here:
U1_p = zeros(JMAX, IMAX);
U2_p = zeros(JMAX, IMAX);
U3_p = zeros(JMAX, IMAX);
U5_p = zeros(JMAX, IMAX);
rho_p = zeros(JMAX, IMAX);
u_p = zeros(JMAX, IMAX);
v_p = zeros(JMAX, IMAX);
T_p = zeros(JMAX, IMAX);
p_p = zeros(JMAX, IMAX);
% Start the main loop
converged = false;
rho_old = rho;
while (~converged && t <= MAXIT)
    % The first thing to do is to compute the appropriate value of the time step needed to satisfy
    % the CFL stability criterion. Only internal points are used in the calculation.
    v_prime = max(max(4/3*mu(2:JMAX-1,2:IMAX-1).^2*gamma./(Pr*rho(2:JMAX-1,2:IMAX-1))));
    delta_t_CFL = 1./(abs(u(2:JMAX-1,2:IMAX-1))/DX + abs(v(2:JMAX-1,2:IMAX-1))/DY + ...
        sqrt(gamma*R*T(2:JMAX-1,2:IMAX-1))*sqrt(1/DX^2 + 1/DY^2) + 2*v_prime*(1/DX^2 + 1/DY^2));
    delta_t(t) = CFL*min(min(delta_t_CFL));

    % Now we are ready to apply MacCormack's technique
    % The first stage is to compute the solution vector U (derivatives wrt t) from the primitive flow field variables
    U1 = rho; % From continuity
    U2 = rho.*u; % From momentum-x
    U3 = rho.*v; % From momentum-y
    U5 = rho.*(c_v*T + (u.^2 + v.^2)/2); % From energy

    %%%%%%%%%%%%%%%%%%
    % Predictor step %
    %%%%%%%%%%%%%%%%%%

    % On first place we have to calculate the flux vectors E (derivatives wrt x) and F (derivatives wrt y) from the primitive variables
    [E1, E2, E3, E5] = Primitive2E(rho, u, p, v, T, mu, lambda, k, c_v, DX, DY, 'Predict_E');
    [F1, F2, F3, F5] = Primitive2F(rho, u, p, v, T, mu, lambda, k, c_v, DX, DY, 'Predict_F');

    % With this we are ready to predict the flow field properties at the next time step

    % For interior points only:
    % Forward differences
    for i = 2:IMAX-1
        for j = 2:JMAX-1
            U1_p(j,i) = U1(j,i) - delta_t(t)*((E1(j,i+1) - E1(j,i))/DX + (F1(j+1,i) - F1(j,i))/DY);
            U2_p(j,i) = U2(j,i) - delta_t(t)*((E2(j,i+1) - E2(j,i))/DX + (F2(j+1,i) - F2(j,i))/DY);
            U3_p(j,i) = U3(j,i) - delta_t(t)*((E3(j,i+1) - E3(j,i))/DX + (F3(j+1,i) - F3(j,i))/DY);
            U5_p(j,i) = U5(j,i) - delta_t(t)*((E5(j,i+1) - E5(j,i))/DX + (F5(j+1,i) - F5(j,i))/DY);
        end
    end

    % Now, decode the flow field variables that will be needed for the calculation of predicted values of the flux vectors E and F.
    % The variables needed are: rho_p, u_p, v_p, p_p, T_p, mu_p, lambda_p and k_p.
    % For interior points we decode them from the U_p values predicted above:
    [rho_p(2:JMAX-1,2:IMAX-1), u_p(2:JMAX-1,2:IMAX-1), v_p(2:JMAX-1,2:IMAX-1), T_p(2:JMAX-1,2:IMAX-1)] = U2Primitive(U1_p(2:JMAX-1,2:IMAX-1), ...
        U2_p(2:JMAX-1,2:IMAX-1), U3_p(2:JMAX-1,2:IMAX-1), U5_p(2:JMAX-1,2:IMAX-1), c_v);
    p_p(2:JMAX-1,2:IMAX-1) = rho_p(2:JMAX-1,2:IMAX-1)*R.*T_p(2:JMAX-1,2:IMAX-1);

    % At the boundary points apply the appropriate boundary conditions:
    [rho_p, u_p, v_p, p_p, T_p] = BC(rho_p, u_p, v_p, p_p, T_p, rho_inf, M_inf*a_inf, p_inf, T_inf, [], R, true);

    % Then, for the remaining properties we have:
    mu_p = DYNVIS(T_p, mu_0, T_0);
    lambda_p = - 2/3*mu_p;
    k_p = THERMC(Pr, c_p, mu_p);

    % This is the end of the predictor step. The necessary flow field properties are now predicted at each grid point.
    % Now we move to the corrector step:

    %%%%%%%%%%%%%%%%%%
    % Corrector step %
    %%%%%%%%%%%%%%%%%%

    % As before, first of all we have to calculate the flux vectors E (derivatives wrt x) and F (derivatives wrt y) from the primitive variables.
    % In this case from the predicted variables we have:
    [E1_p, E2_p, E3_p, E5_p] = Primitive2E(rho_p, u_p, p_p, v_p, T_p, mu_p, lambda_p, k_p, c_v, DX, DY, 'Correct_E');
    [F1_p, F2_p, F3_p, F5_p] = Primitive2F(rho_p, u_p, p_p, v_p, T_p, mu_p, lambda_p, k_p, c_v, DX, DY, 'Correct_F');

    % With this we are ready to correct the flow field properties at the next time step

    % For interior points only:
    % Rearward differences
    for i = 2:IMAX-1
        for j = 2:JMAX-1
            U1(j,i) = 1/2*(U1(j,i) + U1_p(j,i) - delta_t(t)*((E1_p(j,i) - E1_p(j,i-1))/DX + (F1_p(j,i) - F1_p(j-1,i))/DY));
            U2(j,i) = 1/2*(U2(j,i) + U2_p(j,i) - delta_t(t)*((E2_p(j,i) - E2_p(j,i-1))/DX + (F2_p(j,i) - F2_p(j-1,i))/DY));
            U3(j,i) = 1/2*(U3(j,i) + U3_p(j,i) - delta_t(t)*((E3_p(j,i) - E3_p(j,i-1))/DX + (F3_p(j,i) - F3_p(j-1,i))/DY));
            U5(j,i) = 1/2*(U5(j,i) + U5_p(j,i) - delta_t(t)*((E5_p(j,i) - E5_p(j,i-1))/DX + (F5_p(j,i) - F5_p(j-1,i))/DY));
        end
    end

    % Finally, decode the corrected flow field variables.
    % For interior points we decode them from the corrected U values obtained above:
    [rho(2:JMAX-1,2:IMAX-1), u(2:JMAX-1,2:IMAX-1), v(2:JMAX-1,2:IMAX-1), T(2:JMAX-1,2:IMAX-1)] = U2Primitive(U1(2:JMAX-1,2:IMAX-1), ...
        U2(2:JMAX-1,2:IMAX-1), U3(2:JMAX-1,2:IMAX-1), U5(2:JMAX-1,2:IMAX-1), c_v);
    p(2:JMAX-1,2:IMAX-1) = rho(2:JMAX-1,2:IMAX-1)*R.*T(2:JMAX-1,2:IMAX-1);

    % At the boundary points apply the appropriate boundary conditions:
    [rho, u, v, p, T] = BC(rho, u, v, p, T, rho_inf, M_inf*a_inf, p_inf, T_inf, [], R, true);

    % Then, for the remaining properties we have:
    mu = DYNVIS(T, mu_0, T_0);
    lambda = - 2/3*mu;
    k = THERMC(Pr, c_p, mu);

    % This is the end of the corrector step. The necessary flow field properties are now known at each grid point at time t+1.
    % At this point we end the MacCormack's algorithm.

    %%%%%%%%%%%%%%%%%%%%%
    % Convergence check %
    %%%%%%%%%%%%%%%%%%%%%

    % Check for convergence before advancing to the next iteration (include NaN, complex check etc.)
    converged = CONVER(rho_old, rho);
    rho_old = rho;

    time = time + delta_t(t);
    t = t + 1
end
if (converged)
    continuity = MDOT(rho, u, rho_inf, M_inf*a_inf, LVERT, DY); % If convergence has been achieved, check the validity of the solution by confirming conservation of mass
else
    error('CALCULATION FAILED: the maximum number of iterations has been reached before achieving convergence.')
end
if (continuity)
    disp('The calculation for M = 4 with an adiabatic wall boundary condition has finished successfully.')
else
    error('The solution has converged to an invalid result.')
end
M = sqrt(u.^2 + v.^2)./sqrt(gamma*R*T);
figure
surf(x,y,M)
view(2)
axis equal
colorbar
title('Mach number (adiabatic wall)')
xlabel('x')
ylabel('y')
