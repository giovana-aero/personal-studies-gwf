clc,clear

% v2 changes: further modifications to the discretizations (tau terms and such)

% Noteworthy errors:
% - used lambda instead of mu in tau_xy formulas
% - used matrix operations instead of elementwise ones in calc_Et and in calc_delta_t
%   (as for the second, a reminder: ^-1 inverts matrices)
% - forgot minus sign in qx and qy equations
% - wrong values of cp and cv

tic

% Thermodynamics data
M_air = 0.028964; % [kg/mol]
R = 8.314462/M_air; % J/(mol*k)
gamma = 1.4;
Pr = 0.71; % Prandtl number
cv = R/(gamma-1);
cp = gamma*cv;
%cp = 1.003; % (https://www.ohio.edu/mechanical/thermo/property_tables/air/air_Cp_Cv.html)
%cv = 0.716;
%lambda = -2/3*mu; (check page 83)

% Airflow and geometry data (air at sea level conditions)
M1 = 4;
p1 = 1.01325e5; % [Pa]
T1 = 288.16; % [K]
a1 = sqrt(gamma*R*T1);
u1 = a1*M1;
v1 = 0;
rho1 = p1/(R*T1); % [kg/m^3]
mu1 = 1.789e-5; % [kg/m/s] (http://www.aerodynamics4students.com/properties-of-the-atmosphere/sea-level-conditions.php)
L = 0.00001; % Length of the plate [m]
Re = u1*L*rho1/mu1;
H = 5*(5*L/sqrt(Re)); % Height of the computational plane [m]
Tw = T1; % Temperature of the plate at the surface [K]

% Solution configuration
K = 0.6; % Courant number
iter = 4000;
tol = 1e-5; % Tolerance value to check convergence
VIS = 1; % Turn viscosity on and off (probably doesn't work)

% Mesh data
NX = 70;
NY = 70;
delta_x = L/(NX-1);
delta_y = H/(NY-1);

% Initialize solution
u_mat = ones(NY,NX)*u1;
v_mat = zeros(NY,NX);
a_mat = ones(NY,NX)*a1;
T_mat = ones(NY,NX)*T1;
mu_mat = ones(NY,NX)*mu1;
rho_mat = ones(NY,NX)*rho1;
p_mat = ones(NY,NX)*p1;
e_mat = ones(NY,NX)*cv*T1;
%k_mat = ones(NY,NX);
u_bar = zeros(NY,NX);
v_bar = zeros(NY,NX);
a_bar = zeros(NY,NX);
T_bar = zeros(NY,NX);
mu_bar = zeros(NY,NX);
rho_bar = zeros(NY,NX);
p_bar = zeros(NY,NX);
e_bar = zeros(NY,NX);
U1 = ones(NY,NX);
U2 = ones(NY,NX);
U3 = ones(NY,NX);
U4 = ones(NY,NX);
E1 = ones(NY,NX);
E2 = ones(NY,NX);
E3 = ones(NY,NX);
E4 = ones(NY,NX);
F1 = ones(NY,NX);
F2 = ones(NY,NX);
F3 = ones(NY,NX);
F4 = ones(NY,NX);
U1_bar = zeros(NY,NX);
U2_bar = zeros(NY,NX);
U3_bar = zeros(NY,NX);
U4_bar = zeros(NY,NX);
E1_bar = zeros(NY,NX);
E2_bar = zeros(NY,NX);
E3_bar = zeros(NY,NX);
E4_bar = zeros(NY,NX);
F1_bar = zeros(NY,NX);
F2_bar = zeros(NY,NX);
F3_bar = zeros(NY,NX);
F4_bar = zeros(NY,NX);
diff_U1 = zeros(NY,NX);
diff_U2 = zeros(NY,NX);
diff_U3 = zeros(NY,NX);
diff_U4 = zeros(NY,NX);

% Load functions
c10_3_supersonic_flat_plate_functions_v3;

% Boundary conditions
[u_mat,v_mat,T_mat,p_mat,rho_mat,e_mat,mu_mat] = boundaries(u_mat,v_mat,T_mat,p_mat,rho_mat,e_mat,mu_mat,u1,v1,T1,p1,rho1,mu1,cv,R,Tw,NX,NY);

% Calculate other properties from known data
a_mat = sqrt(gamma*R*T_mat);
k_mat = calc_k(Pr,mu_mat,cp);

% Initialize U matrices
[U1,U2,U3,U4] = U_from_air(u_mat,v_mat,rho_mat,e_mat);

% Store residuals
res_his = zeros(4,iter);
res_x = 2;
res_y = 2;

% Check how much time has elapsed in the simulation
sim_time = 0;

% Main loop
for loop = 1:iter

    disp(['Iteration ',num2str(loop)])

    % Keep old values of rho to check convergence
    rho_old = rho_mat;

    % Obtain time step
    delta_t = calc_delta_t(u_mat,v_mat,a_mat,rho_mat,mu_mat,delta_x,delta_y,gamma,Pr,K);
    sim_time = sim_time + delta_t;

    % MacCormack: predictor step
    tau_xx = calc_tau_xx_P(u_mat,v_mat,mu_mat,delta_x,delta_y,NX,NY);
    tau_yy = calc_tau_yy_P(u_mat,v_mat,mu_mat,delta_x,delta_y,NX,NY);
    tau_xy_E = calc_tau_xy_EP(u_mat,v_mat,mu_mat,delta_x,delta_y,NX,NY);
    tau_xy_F = calc_tau_xy_FP(u_mat,v_mat,mu_mat,delta_x,delta_y,NX,NY);
    qx = calc_qx_P(T_mat,k_mat,delta_x,NX,NY);
    qy = calc_qy_P(T_mat,k_mat,delta_y,NX,NY);

    % Calculate relevant values of E and F
    [E1,E2,E3,E4] = E_from_air(u_mat,v_mat,p_mat,rho_mat,e_mat,tau_xx,tau_xy_E,qx);
    [F1,F2,F3,F4] = F_from_air(u_mat,v_mat,p_mat,rho_mat,e_mat,tau_yy,tau_xy_F,qy);

    % Obtain estimated values of U
    for i = 2:NX-1
        for j = 2:NY-1
            U1_bar(j,i) = U1(j,i) - delta_t/delta_x*(E1(j,i+1)-E1(j,i)) - delta_t/delta_y*(F1(j+1,i)-F1(j,i));
            U2_bar(j,i) = U2(j,i) - delta_t/delta_x*(E2(j,i+1)-E2(j,i)) - delta_t/delta_y*(F2(j+1,i)-F2(j,i));
            U3_bar(j,i) = U3(j,i) - delta_t/delta_x*(E3(j,i+1)-E3(j,i)) - delta_t/delta_y*(F3(j+1,i)-F3(j,i));
            U4_bar(j,i) = U4(j,i) - delta_t/delta_x*(E4(j,i+1)-E4(j,i)) - delta_t/delta_y*(F4(j+1,i)-F4(j,i));

            % First step of residual calculations
            temp = [-1/delta_x*(E1(j,i+1)-E1(j,i)) - 1/delta_y*(F1(j+1,i)-F1(j,i));
                    -1/delta_x*(E2(j,i+1)-E2(j,i)) - 1/delta_y*(F2(j+1,i)-F2(j,i));
                    -1/delta_x*(E3(j,i+1)-E3(j,i)) - 1/delta_y*(F3(j+1,i)-F3(j,i));
                    -1/delta_x*(E4(j,i+1)-E4(j,i)) - 1/delta_y*(F4(j+1,i)-F4(j,i))];
            res_his(:,loop) = temp;
        end
    end

    % Obtain estimated air properties
    [u_bar(2:end-1,2:end-1),v_bar(2:end-1,2:end-1),rho_bar(2:end-1,2:end-1),p_bar(2:end-1,2:end-1),T_bar(2:end-1,2:end-1),e_bar(2:end-1,2:end-1)] = air_from_U(U1_bar(2:end-1,2:end-1),U2_bar(2:end-1,2:end-1),U3_bar(2:end-1,2:end-1),U4_bar(2:end-1,2:end-1),cv,R);
    mu_bar(2:end-1,2:end-1) = calc_mu(T_bar(2:end-1,2:end-1),T1,mu1);
    k_bar = calc_k(Pr,mu_bar,cp);

    % Apply boundary conditions to the estimated values
    [u_bar,v_bar,T_bar,p_bar,rho_bar,e_bar,mu_bar] = boundaries(u_bar,v_bar,T_bar,p_bar,rho_bar,e_bar,mu_bar,u1,v1,T1,p1,rho1,mu1,cv,R,Tw,NX,NY);

    % MacCormack: corrector step
    tau_xx = calc_tau_xx_C(u_bar,v_bar,mu_bar,delta_x,delta_y,NX,NY);
    tau_yy = calc_tau_yy_C(u_bar,v_bar,mu_bar,delta_x,delta_y,NX,NY);
    tau_xy_E = calc_tau_xy_EC(u_bar,v_bar,mu_bar,delta_x,delta_y,NX,NY);
    tau_xy_F = calc_tau_xy_FC(u_bar,v_bar,mu_bar,delta_x,delta_y,NX,NY);
    qx = calc_qx_C(T_bar,k_bar,delta_x,NX,NY);
    qy = calc_qy_C(T_bar,k_bar,delta_y,NX,NY);

    % Calculate relevant values of E and F
    [E1_bar,E2_bar,E3_bar,E4_bar] = E_from_air(u_bar,v_bar,p_bar,rho_bar,e_bar,tau_xx,tau_xy_E,qx);
    [F1_bar,F2_bar,F3_bar,F4_bar] = F_from_air(u_bar,v_bar,p_bar,rho_bar,e_bar,tau_yy,tau_xy_F,qy);

    % Obtain corrected values of U
    for i = 2:NX-1
        for j = 2:NY-1
            U1(j,i) = (U1(j,i) + U1_bar(j,i) - delta_t/delta_x*(E1_bar(j,i)-E1_bar(j,i-1)) - delta_t/delta_y*(F1_bar(j,i)-F1_bar(j-1,i)))/2;
            U2(j,i) = (U2(j,i) + U2_bar(j,i) - delta_t/delta_x*(E2_bar(j,i)-E2_bar(j,i-1)) - delta_t/delta_y*(F2_bar(j,i)-F2_bar(j-1,i)))/2;
            U3(j,i) = (U3(j,i) + U3_bar(j,i) - delta_t/delta_x*(E3_bar(j,i)-E3_bar(j,i-1)) - delta_t/delta_y*(F3_bar(j,i)-F3_bar(j-1,i)))/2;
            U4(j,i) = (U4(j,i) + U4_bar(j,i) - delta_t/delta_x*(E4_bar(j,i)-E4_bar(j,i-1)) - delta_t/delta_y*(F4_bar(j,i)-F4_bar(j-1,i)))/2;

             % Second step of residuals calculations
            temp = [-delta_x*(E1_bar(j,i)-E1_bar(j,i-1)) - delta_y*(F1_bar(j,i)-F1_bar(j-1,i));
                    -delta_x*(E2_bar(j,i)-E2_bar(j,i-1)) - delta_y*(F1_bar(j,i)-F1_bar(j-1,i));
                    -delta_x*(E3_bar(j,i)-E3_bar(j,i-1)) - delta_y*(F1_bar(j,i)-F1_bar(j-1,i));
                    -delta_x*(E4_bar(j,i)-E4_bar(j,i-1)) - delta_y*(F1_bar(j,i)-F1_bar(j-1,i))];
            res_his(:,loop) = (res_his(:,loop) + temp)/2;
        end
    end

    % Calculate air properties at the interior
    [u_mat(2:end-1,2:end-1),v_mat(2:end-1,2:end-1),rho_mat(2:end-1,2:end-1),p_mat(2:end-1,2:end-1),T_mat(2:end-1,2:end-1),e_mat(2:end-1,2:end-1)] = air_from_U(U1(2:end-1,2:end-1),U2(2:end-1,2:end-1),U3(2:end-1,2:end-1),U4(2:end-1,2:end-1),cv,R);
    mu_mat = calc_mu(T_mat,T1,mu1);

    % Apply boundary conditions
    [u_mat,v_mat,T_mat,p_mat,rho_mat,e_mat,mu_mat] = boundaries(u_mat,v_mat,T_mat,p_mat,rho_mat,e_mat,mu_mat,u1,v1,T1,p1,rho1,mu1,cv,R,Tw,NX,NY);
    a_mat = sqrt(gamma*R*T_mat);

    % Obtain U values at borders 3 and 4
    [U1(1,2:end),U2(1,2:end),U3(1,2:end),U4(1,2:end)] = U_from_air(u_mat(1,2:end),v_mat(1,2:end),rho_mat(1,2:end),e_mat(1,2:end)); % Border 3
    [U1(2:end-1,end),U2(2:end-1,end),U3(2:end-1,end),U4(2:end-1,end)] = U_from_air(u_mat(2:end-1,end),v_mat(2:end-1,end),rho_mat(2:end-1,end),e_mat(2:end-1,end)); % Border 4

    % Obtain U values in the interior
    [U1(2:end-1,2:end-1),U2(2:end-1,2:end-1),U3(2:end-1,2:end-1),U4(2:end-1,2:end-1)] = U_from_air(u_mat(2:end-1,2:end-1),v_mat(2:end-1,2:end-1),rho_mat(2:end-1,2:end-1),e_mat(2:end-1,2:end-1));

    % Check convergence
    status = convergence(rho_mat,rho_old,tol);
    if status == 1
        disp('Simulation convergence')
        break
    end

    if ~isreal(mu_mat)
        error('Broken solution - presence of complex numbers')
    end
end

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
t = toc;
min = fix(t/60); s = rem(t,60);
fprintf('Time: %d min, %.2f s\n', min, s)
%min = fix(sim_time/60); s = rem(sim_time,60);
%fprintf('Simulation time: %d min ,%f s\n', min, s)

% Plot contours
figure(1),clf
M_mat = sqrt(u_mat.^2+v_mat.^2)./a_mat;
contourf(M_mat)
title('Mach number')
colorbar
colormap turbo

figure(2),clf
contourf(T_mat)
title('Temperature [K]')
colorbar
colormap jet

figure(3),clf
contourf(p_mat)
title('Pressure [Pa]')
colorbar
colormap jet

figure(4),clf
contourf(rho_mat)
title('Density [kg/m^3]')
colorbar
colormap jet

% Plot residuals
res_his = abs(res_his);
figure(5),clf
subplot(2,2,1)
plot(1:loop,res_his(1,1:loop)),grid on
xlabel('Iterations')
ylabel('Residuals')
title('U1')
subplot(2,2,2)
plot(1:loop,res_his(2,1:loop)),grid on
xlabel('Iterations')
ylabel('Residuals')
title('U2')
subplot(2,2,3)
plot(1:loop,res_his(3,1:loop)),grid on
xlabel('Iterations')
ylabel('Residuals')
title('U3')
subplot(2,2,4)
plot(1:loop,res_his(4,1:loop)),grid on
xlabel('Iterations')
ylabel('Residuals')
title('U4')


%MDOT(rho_mat, u_mat, rho1, u1, H, delta_y)
