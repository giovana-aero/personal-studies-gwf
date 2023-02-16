clc,clear

% v2 changes: modifications to the discretizations (tau terms and such)

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

% Solution configuration
K = 0.6; % Courant number
iter = 5;
VIS = 1; % Turn viscosity on and off (probably doesn't work)

% Mesh data
NX = 60;
NY = 60;
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
c10_3_supersonic_flat_plate_functions;

% Boundary conditions
[u_mat,v_mat,T_mat,p_mat,rho_mat,e_mat,mu_mat] = boundaries(u_mat,v_mat,T_mat,p_mat,rho_mat,e_mat,mu_mat,u1,v1,T1,p1,rho1,mu1,cv,R,Tw,NX,NY);


%rho_mat(1:end-1,2:end) = rho_mat(1:end-1,2:end) + rho1; is this necessary?



% Calculate other properties from known data
%mu_mat = calc_mu(T_mat,T1,mu1);
%e_mat = cv*T_mat;
a_mat = sqrt(gamma*R*T_mat);

% Initialize U matrices
[U1,U2,U3,U4] = U_from_air(u_mat,v_mat,rho_mat,e_mat);


%rho_mat(:,1) = p_mat(:,1)./(R*T_mat(:,1));
%rho_mat(end,:) = rho_mat(end,:) + rho1;
%k_mat = calc_k(Pr,mu_mat,cp);

% Main loop
for loop = 1:iter

    disp(['Iteration ',num2str(loop)])

    % Obtain time step
    delta_t = calc_delta_t(u_mat,v_mat,a_mat,rho_mat,mu_mat,delta_x,delta_y,gamma,Pr,K);

    % MacCormack: predictor step
    for i = 2:NX-1
        for j = 2:NY-1
%                if j == 5 && i == 2
%                    1;
%                    break
%                end
            % Obtain relevant terms to use in En_ipj
            lambda = -2/3*mu_mat(j,i+1);
            tau_xx = lambda*((u_mat(j,i+1)-u_mat(j,i))/delta_x+(v_mat(j+1,i+1)-v_mat(j-1,i+1))/(2*delta_y)) + 2*mu_mat(j,i+1)*(u_mat(j,i+1)-u_mat(j,i))/delta_x*VIS;
            tau_xy = mu_mat(j,i+1)*((u_mat(j+1,i+1)-u_mat(j-1,i+1))/(2*delta_y)+(v_mat(j,i+1)-v_mat(j,i))/delta_x)*VIS;
            k = calc_k(Pr,mu_mat(j,i+1),cp);
            qx = -k*(T_mat(j,i+1)-T_mat(j,i))/delta_x;
            Et = calc_Et(u_mat(j,i+1),v_mat(j,i+1),rho_mat(j,i+1),e_mat(j,i+1));

            % Calculate E1_ipj
            E1_ipj = rho_mat(j,i+1)*u_mat(j,i+1);
            % Calculate E2_ipj
            E2_ipj = rho_mat(j,i+1)*u_mat(j,i+1)^2 + p_mat(j,i+1) - tau_xx;
            % Calculate E3_ipj
            E3_ipj = rho_mat(j,i+1)*u_mat(j,i+1)*v_mat(j,i+1) - tau_xy;
            % Calculate E4_ipj
            E4_ipj = (Et+p_mat(j,i+1))*u_mat(j,i+1) - u_mat(j,i+1)*tau_xx - v_mat(j,i+1)*tau_xy + qx;

            % Obtain relevant terms to use in En_ij
            lambda = -2/3*mu_mat(j,i);
            tau_xx = lambda*((u_mat(j,i)-u_mat(j,i-1))/delta_x+(v_mat(j+1,i)-v_mat(j-1,i))/(2*delta_y)) + 2*mu_mat(j,i)*(u_mat(j,i)-u_mat(j,i-1))/delta_x*VIS;
            tau_xy = mu_mat(j,i)*((u_mat(j+1,i)-u_mat(j-1,i))/(2*delta_y)+(v_mat(j,i)-v_mat(j,i-1))/delta_x)*VIS;
            k = calc_k(Pr,mu_mat(j,i),cp);
            qx = -k*(T_mat(j,i)-T_mat(j,i-1))/delta_x;
            Et = calc_Et(u_mat(j,i),v_mat(j,i),rho_mat(j,i),e_mat(j,i));

            % Calculate E1_ij
            E1_ij = rho_mat(j,i)*u_mat(j,i);
            % Calculate E2_ij
            E2_ij = rho_mat(j,i)*u_mat(j,i)^2 + p_mat(j,i) - tau_xx;
            % Calculate E3_ij
            E3_ij = rho_mat(j,i)*u_mat(j,i)*v_mat(j,i) - tau_xy;
            % Calculate E4_ij
            E4_ij = (Et+p_mat(j,i))*u_mat(j,i) - u_mat(j,i)*tau_xx - v_mat(j,i)*tau_xy + qx;

            % Obtain relevant terms to use in Fn_ijp
            lambda = -2/3*mu_mat(j+1,i);
            tau_yy = lambda*((u_mat(j+1,i+1)-u_mat(j+1,i-1))/(2*delta_x)+(v_mat(j+1,i)-v_mat(j,i))/delta_y) + 2*mu_mat(j+1,i)*(v_mat(j+1,i)-v_mat(j,i))/delta_y*VIS;
            tau_xy = mu_mat(j+1,i)*((u_mat(j+1,i)-u_mat(j,i))/delta_y+(v_mat(j+1,i+1)-v_mat(j+1,i-1))/(2*delta_x))*VIS;
            k = calc_k(Pr,mu_mat(j+1,i),cp);
            qy = -k*(T_mat(j+1,i)-T_mat(j,i))/delta_y;
            Et = calc_Et(u_mat(j+1,i),v_mat(j+1,i),rho_mat(j+1,i),e_mat(j+1,i));

            % Calculate F1_ijp
            F1_ijp = rho_mat(j+1,i)*v_mat(j+1,i);
            % Calculate F2_ijp
            F2_ijp = rho_mat(j+1,i)*u_mat(j+1,i)*v_mat(j+1,i) - tau_xy;
            % Calculate F3_ijp
            F3_ijp = rho_mat(j+1,i)*v_mat(j+1,i)^2 + p_mat(j+1,i) - tau_yy;
            % Calculate F4_ijp
            F4_ijp = (Et+p_mat(j+1,i))*v_mat(j+1,i) - u_mat(j+1,i)*tau_xy - v_mat(j+1,i)*tau_yy + qy;

            % Obtain relevant terms to use in Fn_ij
            lambda = -2/3*mu_mat(j,i);
            tau_yy = lambda*((u_mat(j,i+1)-u_mat(j,i-1))/(2*delta_x)+(v_mat(j,i)-v_mat(j-1,i))/delta_y) + 2*mu_mat(j,i)*(v_mat(j,i)-v_mat(j-1,i))/delta_y*VIS;
            tau_xy = mu_mat(j,i)*((u_mat(j,i)-u_mat(j-1,i))/delta_y+(v_mat(j,i+1)-v_mat(j,i-1))/(2*delta_x))*VIS;
            k = calc_k(Pr,mu_mat(j,i),cp);
            qy = -k*(T_mat(j,i)-T_mat(j-1,i))/delta_y;
            Et = calc_Et(u_mat(j,i),v_mat(j,i),rho_mat(j,i),e_mat(j,i));

            % Calculate F1_ij
            F1_ij = rho_mat(j,i)*v_mat(j,i);
            % Calculate F2_ij
            F2_ij = rho_mat(j,i)*u_mat(j,i)*v_mat(j,i) - tau_xy;
            % Calculate F3_ij
            F3_ij = rho_mat(j,i)*v_mat(j,i)^2 + p_mat(j,i) - tau_yy;
            % Calculate F4_ij
            F4_ij = (Et+p_mat(j,i))*v_mat(j,i) - u_mat(j,i)*tau_xy - v_mat(j,i)*tau_yy + qy;

            % Estimate U terms
            U1_bar(j,i) = U1(j,i) - delta_t/delta_x*(E1_ipj-E1_ij) - delta_t/delta_y*(F1_ijp-F1_ij);
            U2_bar(j,i) = U2(j,i) - delta_t/delta_x*(E2_ipj-E2_ij) - delta_t/delta_y*(F2_ijp-F2_ij);
            U3_bar(j,i) = U3(j,i) - delta_t/delta_x*(E3_ipj-E3_ij) - delta_t/delta_y*(F3_ijp-F3_ij);
            U4_bar(j,i) = U4(j,i) - delta_t/delta_x*(E4_ipj-E4_ij) - delta_t/delta_y*(F4_ijp-F4_ij);

        end
    end

    % Obtain estimated air properties
    [u_bar(2:end-1,2:end-1),v_bar(2:end-1,2:end-1),rho_bar(2:end-1,2:end-1),p_bar(2:end-1,2:end-1),T_bar(2:end-1,2:end-1),e_bar(2:end-1,2:end-1)] = air_from_U(U1_bar(2:end-1,2:end-1),U2_bar(2:end-1,2:end-1),U3_bar(2:end-1,2:end-1),U4_bar(2:end-1,2:end-1),cv,R);
    mu_bar(2:end-1,2:end-1) = calc_mu(T_bar(2:end-1,2:end-1),T1,mu1);

%    % Remove small values of temperature (remove later?)
%    for i = 2:NX-1
%        for j = 2:NY-1
%            if abs(T_bar(j,i)) <= 1e-10
%                T_bar(j,i) = 0;
%             end
%        end
%    end

    % Apply boundary conditions to the estimated values
    [u_bar,v_bar,T_bar,p_bar,rho_bar,e_bar,mu_bar] = boundaries(u_bar,v_bar,T_bar,p_bar,rho_bar,e_bar,mu_bar,u1,v1,T1,p1,rho1,mu1,cv,R,Tw,NX,NY);



%    mu_bar = real(mu_bar); % remove later?



    % MacCormack: corrector step
    for i = 2:NX-1
        for j = 2:NY-1

            % Obtain relevant terms to use in En_ij
            lambda = -2/3*mu_bar(j,i);
            tau_xx = lambda*((u_bar(j,i+1)-u_bar(j,i))/delta_x+(v_bar(j+1,i)-v_bar(j-1,i))/(2*delta_y)) + 2*mu_bar(j,i)*(u_bar(j,i+1)-u_bar(j,i))/delta_x;
            tau_xy = mu_bar(j,i)*((u_bar(j+1,i)-u_bar(j-1,i))/(2*delta_y)+(v_bar(j,i+1)-v_bar(j,i))/delta_x);
            k = calc_k(Pr,mu_bar(j,i),cp);
            qx = -k*(T_bar(j,i+1)-T_bar(j,i))/delta_x;
            Et = calc_Et(u_bar(j,i),v_bar(j,i),rho_bar(j,i),e_bar(j,i));

            % Calculate E1_ij
            E1_bar_ij = rho_bar(j,i)*u_bar(j,i);
            % Calculate E2_ij
            E2_bar_ij = rho_bar(j,i)*u_bar(j,i)^2 + p_bar(j,i) - tau_xx;
            % Calculate E3_ij
            E3_bar_ij = rho_bar(j,i)*u_bar(j,i)*v_bar(j,i) - tau_xy;
            % Calculate E4_ij
            E4_bar_ij = (Et+p_bar(j,i))*u_bar(j,i) - u_bar(j,i)*tau_xx - v_bar(j,i)*tau_xy + qx;

            % Obtain relevant terms to use in En_imj
            lambda = -2/3*mu_bar(j,i-1);
            tau_xx = lambda*((u_bar(j,i)-u_bar(j,i-1))/delta_x+(v_bar(j+1,i-1)-v_bar(j-1,i-1))/(2*delta_y)) + 2*mu_bar(j,i-1)*(u_bar(j,i)-u_bar(j,i-1))/delta_x;
            tau_xy = mu_bar(j,i-1)*((u_bar(j+1,i-1)-u_bar(j-1,i-1))/(2*delta_y)+(v_bar(j,i)-v_bar(j,i-1))/delta_x);
            k = calc_k(Pr,mu_bar(j,i-1),cp);
            qx = -k*(T_bar(j,i)-T_bar(j,i-1))/delta_x;
            Et = calc_Et(u_bar(j,i-1),v_bar(j,i-1),rho_bar(j,i-1),e_bar(j,i-1));

            % Calculate E1_imj
            E1_bar_imj = rho_bar(j,i-1)*u_bar(j,i-1);
            % Calculate E2_imj
            E2_bar_imj = rho_bar(j,i-1)*u_bar(j,i-1)^2 + p_bar(j,i-1) - tau_xx;
            % Calculate E3_imj
            E3_bar_imj = rho_bar(j,i-1)*u_bar(j,i-1)*v_bar(j,i-1) - tau_xy;
            % Calculate E4_imj
            E4_bar_imj = (Et+p_bar(j,i-1))*u_bar(j,i-1) - u_bar(j,i-1)*tau_xx - v_bar(j,i-1)*tau_xy + qx;

            % Obtain relevant terms to use in Fn_ij
            lambda = -2/3*mu_bar(j,i);
            tau_yy = lambda*((u_bar(j,i+1)-u_bar(j,i-1))/(2*delta_x)+(v_bar(j+1,i)-v_bar(j,i))/delta_y) + 2*mu_bar(j,i)*(v_bar(j+1,i)-v_bar(j,i))/delta_y;
            tau_xy = mu_bar(j,i)*((u_bar(j+1,i)-u_bar(j,i))/delta_y+(v_bar(j,i+1)-v_bar(j,i-1))/(2*delta_x));
            k = calc_k(Pr,mu_bar(j,i),cp);
            qy = -k*(T_bar(j+1,i)-T_bar(j,i))/delta_y;
            Et = calc_Et(u_bar(j,i),v_bar(j,i),rho_bar(j,i),e_bar(j,i));

            % Calculate F1_ij
            F1_bar_ij = rho_bar(j,i)*v_bar(j,i);
            % Calculate F2_ij
            F2_bar_ij = rho_bar(j,i)*u_bar(j,i)*v_bar(j,i) - tau_xy;
            % Calculate F3_ij
            F3_bar_ij = rho_bar(j,i)*v_bar(j,i)^2 + p_bar(j,i) - tau_yy;
            % Calculate F4_ij
            F4_bar_ij = (Et+p_bar(j,i))*v_bar(j,i) - u_bar(j,i)*tau_xy - v_bar(j,i)*tau_yy + qy;

            % Obtain relevant terms to use in Fn_ijm  *****
            lambda = -2/3*mu_bar(j-1,i);
            tau_yy = lambda*((u_bar(j-1,i+1)-u_bar(j-1,i-1))/(2*delta_x)+(v_bar(j,i)-v_bar(j-1,i))/delta_y) + 2*mu_bar(j-1,i)*(v_bar(j,i)-v_bar(j-1,i))/delta_y;
            tau_xy = mu_bar(j-1,i)*((u_bar(j,i)-u_bar(j-1,i))/delta_y+(v_bar(j-1,i+1)-v_bar(j-1,i-1))/(2*delta_x));
            k = calc_k(Pr,mu_bar(j-1,i),cp);
            qy = -k*(T_bar(j,i)-T_bar(j-1,i))/delta_y;
            Et = calc_Et(u_bar(j-1,i),v_bar(j-1,i),rho_bar(j-1,i),e_bar(j-1,i));

            % Calculate F1_ijm
            F1_bar_ijm = rho_bar(j-1,i)*v_bar(j-1,i);
            % Calculate F2_ijm
            F2_bar_ijm = rho_bar(j-1,i)*u_bar(j-1,i)*v_bar(j-1,i) - tau_xy;
            % Calculate F3_ijm
            F3_bar_ijm = rho_bar(j-1,i)*v_bar(j-1,i)^2 + p_bar(j-1,i) - tau_yy;
            % Calculate F4_ijm
            F4_bar_ijm = (Et+p_bar(j-1,i))*v_bar(j-1,i) - u_bar(j-1,i)*tau_xy - v_bar(j-1,i)*tau_yy + qy;

            % Calculate corrected U terms
            U1(j,i) = (U1(j,i) + U1_bar(j,i) - delta_t/delta_x*(E1_bar_ij-E1_bar_imj) - delta_t/delta_y*(F1_bar_ij-F1_bar_ijm))/2;
            U2(j,i) = (U2(j,i) + U2_bar(j,i) - delta_t/delta_x*(E2_bar_ij-E2_bar_imj) - delta_t/delta_y*(F2_bar_ij-F2_bar_ijm))/2;
            U3(j,i) = (U3(j,i) + U3_bar(j,i) - delta_t/delta_x*(E3_bar_ij-E3_bar_imj) - delta_t/delta_y*(F3_bar_ij-F3_bar_ijm))/2;
            U4(j,i) = (U4(j,i) + U4_bar(j,i) - delta_t/delta_x*(E4_bar_ij-E4_bar_imj) - delta_t/delta_y*(F4_bar_ij-F4_bar_ijm))/2;
        end
    end

    % Calculate air properties at the interior
    [u_mat(2:end-1,2:end-1),v_mat(2:end-1,2:end-1),rho_mat(2:end-1,2:end-1),p_mat(2:end-1,2:end-1),T_mat(2:end-1,2:end-1),e_mat(2:end-1,2:end-1)] = air_from_U(U1(2:end-1,2:end-1),U2(2:end-1,2:end-1),U3(2:end-1,2:end-1),U4(2:end-1,2:end-1),cv,R);
    mu_mat = calc_mu(T_mat,T1,mu1);

    % Apply boundary conditions
    [u_mat,v_mat,T_mat,p_mat,rho_mat,e_mat,mu_mat] = boundaries(u_mat,v_mat,T_mat,p_mat,rho_mat,e_mat,mu_mat,u1,v1,T1,p1,rho1,mu1,cv,R,Tw,NX,NY);

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

