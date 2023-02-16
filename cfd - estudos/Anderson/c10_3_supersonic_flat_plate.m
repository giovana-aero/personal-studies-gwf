clc,clear

% Thermodynamics data
M_air = 0.028964; % [kg/mol]
R = 8.314462/M_air; % J/(mol*k)
gamma = 1.4;
Pr = 0.71; % Prandtl number
cp = 1.003; % (https://www.ohio.edu/mechanical/thermo/property_tables/air/air_Cp_Cv.html)
cv = 0.716;
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

% Solution data
K = 0.5; % Courant number
iter = 2;

% Mesh data
NX = 6;
NY = 6;
delta_x = L/(NX-1);
delta_y = H/(NY-1);

% Initialize solution
u_mat = ones(NY,NX);
v_mat = zeros(NY,NX);
a_mat = ones(NY,NX);
T_mat = ones(NY,NX);
mu_mat = ones(NY,NX);
rho_mat = ones(NY,NX);
p_mat = ones(NY,NX);
%e_mat = ones(NY,NX);
%k_mat = ones(NY,NX);
u_bar = ones(NY,NX);
v_bar = ones(NY,NX);
a_bar = ones(NY,NX);
T_bar = ones(NY,NX);
mu_bar = ones(NY,NX);
rho_bar = ones(NY,NX);
p_bar = ones(NY,NX);
e_bar = ones(NY,NX);
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
c10_3_supersonic_flat_plate_functions;

% Boundary conditions
% Case 1
u_mat(1,1) = 0;
v_mat(1,1) = 0;
T_mat(1,1) = T1;
p_mat(1,1) = p1;
rho_mat(1,1) = rho1;
% Case 2
u_mat(2:end,1) = u_mat(2:end,1) + u1; % Left boundary
v_mat(2:end,1) = v_mat(2:end,1)-v_mat(2:end,1); % Left boundary
v_mat(end,2:end) = v_mat(end,2:end)-v_mat(end,2:end); % Top boundary
T_mat(2:end,1) = T_mat(2:end,1) + T1; % Left boundary
T_mat(end,2:end) = T_mat(end,2:end) + T1; % Top boundary
p_mat(2:end,1) = p_mat(2:end,1) + p1; % Left boundary
p_mat(end,2:end) = p_mat(end,2:end) + p1; % Top boundary
rho_mat(2:end,1) = rho_mat(2:end,1) + rho1; % Left boundary
rho_mat(end,2:end) = rho_mat(end,2:end) + rho1; % Top boundary
% Case 3
T_mat(1,2:end) = T_mat(1,2:end) + Tw;
u_mat(1,2:end) = u_mat(1,2:end) - u_mat(1,2:end);
v_mat(1,2:end) = v_mat(1,2:end) - v_mat(1,2:end);



rho_mat(1:end-1,2:end) = rho_mat(1:end-1,2:end) + rho1;

%disp(u_mat)
%disp(v_mat)
%disp(T_mat)
%disp(p_mat)
%disp(rho_mat)

% Calculate other properties from known data
mu_mat = calc_mu(T_mat,T1,mu1);
e_mat = cv*T_mat;
a_mat = sqrt(gamma*R*T_mat);

%rho_mat(:,1) = p_mat(:,1)./(R*T_mat(:,1));
%rho_mat(end,:) = rho_mat(end,:) + rho1;
%k_mat = calc_k(Pr,mu_mat,cp);


% Main loop
for loop = 1:iter

    disp(['Iteration ',num2str(loop)])

    % Obtain time step
    delta_t = calc_delta_t(u_mat,v_mat,a_mat,rho_mat,mu_mat,delta_x,delta_y,gamma,Pr,K);

    % Obtain values of E and F along the borders
%    [E1,E2,E3,E4,F1,F2,F3,F4] = E_F_borders(u_mat,v_mat,p_mat,rho_mat,T_mat,e_mat,mu_mat,cp,Pr,delta_x,delta_y,NX,NY);

    % Obtain values of E and F inside the borders
    for j = 2:NY-1
        for i = 2:NX-1
            % Calculate tau terms
            tau_xx = calc_tau_xx(u_mat(j,i),u_mat(j,i-1),v_mat(j+1,i),v_mat(j-1,i),mu_mat(j,i),delta_x,delta_y);
            tau_yy = calc_tau_yy(u_mat(j,i+1),u_mat(j,i-1),v_mat(j,i),v_mat(j-1,i),mu_mat(j,i),delta_x,delta_y);
            tau_xy_E = calc_tau_xy_E(u_mat(j+1,i),u_mat(j-1,i),v_mat(j,i),v_mat(j,i-1),mu_mat(j,i),delta_x,delta_y);
            tau_xy_F = calc_tau_xy_F(u_mat(j,i),u_mat(j-1,i),v_mat(j,i+1),v_mat(j,i-1),mu_mat(j,i),delta_x,delta_y);

            % Calculate qx and qy
%            mu = calc_mu(T_mat(j,i),T1,mu1);
            k = calc_k(Pr,mu_mat(j,i),cp);
            [qx,qy] = calc_q(k,T_mat(j,i),T_mat(j,i-1),T_mat(j,i),T_mat(j-1,i),delta_x,delta_y);

            % Calculate E and F terms
            Et = calc_Et(u_mat(j,i),v_mat(j,i),rho_mat(j,i),e_mat(j,i));
            [E1(j,i),E2(j,i),E3(j,i),E4(j,i)] = calc_E(u_mat(j,i),v_mat(j,i),p_mat(j,i),rho_mat(j,i),tau_xx,tau_xy_E,Et,qx);
            [F1(j,i),F2(j,i),F3(j,i),F4(j,i)] = calc_F(u_mat(j,i),v_mat(j,i),p_mat(j,i),rho_mat(j,i),tau_yy,tau_xy_F,Et,qy);
        end
    end

    % MacCormack: calculate estimated values in the internal points and along
    % the left and bottom borders
    for j = 1:NY-1
        for i = 1:NX-1
            % Predictor step
            diff_U1(j,i) = -(E1(j,i+1)-E1(j,i))/delta_x - (F1(j+1,i)-F1(j,i))/delta_y;
            diff_U2(j,i) = -(E2(j,i+1)-E2(j,i))/delta_x - (F2(j+1,i)-F2(j,i))/delta_y;
            diff_U3(j,i) = -(E3(j,i+1)-E3(j,i))/delta_x - (F3(j+1,i)-F3(j,i))/delta_y;
            diff_U4(j,i) = -(E4(j,i+1)-E4(j,i))/delta_x - (F4(j+1,i)-F4(j,i))/delta_y;

            % Get estimated values of U
            U1_bar(j,i) = U1(j,i) + diff_U1(j,i)*delta_t;
            U2_bar(j,i) = U2(j,i) + diff_U2(j,i)*delta_t;
            U3_bar(j,i) = U3(j,i) + diff_U3(j,i)*delta_t;
            U4_bar(j,i) = U4(j,i) + diff_U4(j,i)*delta_t;

            % Calculate estimated air properties
            [u_bar(j,i),v_bar(j,i),rho_bar(j,i),p_bar(j,i),T_bar(j,i),e_bar(j,i)] = air_from_U(U1_bar(j,i),U2_bar(j,i),U3_bar(j,i),U4_bar(j,i),cv,R);
        end
    end

    % Obtain estimated values of E and F along the borders
%    [E1_bar,E2_bar,E3_bar,E4_bar,F1_bar,F2_bar,F3_bar,F4_bar] = E_F_borders(u_bar,v_bar,p_bar,rho_bar,T_bar,e_bar,mu_bar,cp,Pr,delta_x,delta_y,NX,NY);

    % MacCormack: obtain estimated values of E and F inside the borders
    for j = 2:NY-1
        for i = 2:NX-1
            % Calculate tau terms (from estimated data)
            tau_xx = calc_tau_xx(u_bar(j,i+1),u_bar(j,i),v_bar(j+1,i),v_bar(j-1,i),mu_bar(j,i),delta_x,delta_y);
            tau_yy = calc_tau_yy(u_bar(j,i+1),u_bar(j,i-1),v_bar(j+1,i),v_bar(j,i),mu_bar(j,i),delta_x,delta_y);
            tau_xy_E = calc_tau_xy_E(u_bar(j+1,i),u_bar(j-1,i),v_bar(j,i+1),v_bar(j,i),mu_bar(j,i),delta_x,delta_y);
            tau_xy_F = calc_tau_xy_F(u_bar(j+1,i),u_bar(j,i),v_bar(j,i+1),v_bar(j,i-1),mu_bar(j,i),delta_x,delta_y);

            % Calculate qx and qy (from estimated data)
%            mu = calc_mu(T_mat(j,i),T1,mu1);
            k = calc_k(Pr,mu_mat(j,i),cp);
            [qx,qy] = calc_q(k,T_bar(j,i+1),T_bar(j,i),T_bar(j+1,i),T_bar(j,i),delta_x,delta_y);

            % Get estimated values of E and F
            Et = calc_Et(u_bar(j,i),v_bar(j,i),rho_bar(j,i),e_bar(j,i));
            [E1_bar(j,i),E2_bar(j,i),E3_bar(j,i),E4_bar(j,i)] = calc_E(u_bar(j,i),v_bar(j,i),p_bar(j,i),rho_bar(j,i),tau_xx,tau_xy_E,Et,qx);
            [F1_bar(j,i),F2_bar(j,i),F3_bar(j,i),F4_bar(j,i)] = calc_F(u_bar(j,i),v_bar(j,i),p_bar(j,i),rho_bar(j,i),tau_yy,tau_xy_F,Et,qy);
        end
    end

    for j = 2:NY-1
        for i = 2:NX-1
            % Corrector step
            diff_U1_bar = -(E1_bar(j,i)-E1_bar(j,i-1))/delta_x - (F1_bar(j,i)-F1_bar(j-1,i))/delta_y;
            diff_U2_bar = -(E2_bar(j,i)-E2_bar(j,i-1))/delta_x - (F2_bar(j,i)-F2_bar(j-1,i))/delta_y;
            diff_U3_bar = -(E3_bar(j,i)-E3_bar(j,i-1))/delta_x - (F3_bar(j,i)-F3_bar(j-1,i))/delta_y;
            diff_U4_bar = -(E4_bar(j,i)-E4_bar(j,i-1))/delta_x - (F4_bar(j,i)-F4_bar(j-1,i))/delta_y;

            % Average values
            diff_U1_AV = (diff_U1(j,i)+diff_U1_bar)/2;
            diff_U2_AV = (diff_U2(j,i)+diff_U2_bar)/2;
            diff_U3_AV = (diff_U3(j,i)+diff_U3_bar)/2;
            diff_U4_AV = (diff_U4(j,i)+diff_U4_bar)/2;

            % Corrected values
            U1(j,i) = U1(j,i) + diff_U1_AV*delta_t;
            U2(j,i) = U2(j,i) + diff_U2_AV*delta_t;
            U3(j,i) = U3(j,i) + diff_U3_AV*delta_t;
            U4(j,i) = U4(j,i) + diff_U4_AV*delta_t;
        end
    end

    % Calculate air properties at the interior
    [u_mat(2:end-1,2:end-1),v_mat(2:end-1,2:end-1),rho_mat(2:end-1,2:end-1),p_mat(2:end-1,2:end-1),T_mat(2:end-1,2:end-1),e_mat(2:end-1,2:end-1)] = air_from_U(U1(2:end-1,2:end-1),U2(2:end-1,2:end-1),U3(2:end-1,2:end-1),U4(2:end-1,2:end-1),cv,R);

    % Apply boundary conditions
    % Case 3
    p_mat(1,2:end) = 2*p_mat(2,2:end) - p_mat(3,2:end);
    rho_mat(1,2:end) = p_mat(1,2:end)./(R*T_mat(1,2:end));
    % Case 4
    u_mat(2:end-1,NX) = 2*u_mat(2:end-1,NX-1) - u_mat(2:end-1,NX-2);
    v_mat(2:end-1,NX) = 2*v_mat(2:end-1,NX-1) - v_mat(2:end-1,NX-2);
    T_mat(2:end-1,NX) = 2*T_mat(2:end-1,NX-1) - T_mat(2:end-1,NX-2);
    p_mat(2:end-1,NX) = 2*p_mat(2:end-1,NX-1) - p_mat(2:end-1,NX-2);
    rho_mat(2:end-1,NX) = p_mat(2:end-1,NX)./(R*T_mat(2:end-1,NX));
    e_mat(2:end-1,NX) = cv*T_mat(2:end-1,NX);
    mu_mat(2:end-1,NX) = calc_mu(T_mat(2:end-1,NX),T1,mu1);

    % Obtain U values at borders 3 and 4
    [U1(1,2:end),U2(1,2:end),U3(1,2:end),U4(1,2:end)] = U_from_air(u_mat(1,2:end),v_mat(1,2:end),rho_mat(1,2:end),e_mat(1,2:end)); % Border 3
    [U1(2:end-1,end),U2(2:end-1,end),U3(2:end-1,end),U4(2:end-1,end)] = U_from_air(u_mat(2:end-1,end),v_mat(2:end-1,end),rho_mat(2:end-1,end),e_mat(2:end-1,end)); % Border 4

    % Obtain U values in the interior
    [U1(2:end-1,2:end-1),U2(2:end-1,2:end-1),U3(2:end-1,2:end-1),U4(2:end-1,2:end-1)] = U_from_air(u_mat(2:end-1,2:end-1),v_mat(2:end-1,2:end-1),rho_mat(2:end-1,2:end-1),e_mat(2:end-1,2:end-1));

end

disp(U1)
disp(U2)
disp(U3)
disp(U4)

% Obtain final distirbution of sound speed
%a_mat =

