clc,clear

% Thermodynamics data
R = 8.314462; % J/(mol*k)
gamma = 1.4;
M_air = 0.028964; % [kg/mol]

% Flow properties (inlet)
M1 = 2;
p1 = 1.01e5; % [N/m^2]
rho1 = 1.23; % [kg/ m^3]
T1 = 286.1; % [K]
a1 = sqrt(gamma*R*T1/M_air);
u1 = a1*M1;
v1 = 0;

% Geometry data
L = 65; % [m]
E = 10; % [m]
H = 40; % [m]
theta = 5.352; % [degrees]
NY = 41;
delta_eta = 1/(NY-1);
x = 0;

% Solution setup
iter = 100;
C = 0.5;
Cx = 0;

% Calculate height
function h = calc_h(x,E,H,theta)
    if x <= E
        h = H;
    else
        h = H + (x-E)*tand(theta);
    end
end
% Calculate ys
function ys = calc_ys(x,E,theta)
    if x <= E
        ys = 0;
    else
        ys = -(x-E)*tand(theta);
    end
end
% Calculate (d eta)/dx)
function df_eta_x = calc_df_eta_x(x,E,eta,h,theta)
    if x <= E
        df_eta_x = 0;
    else
        df_eta_x = (1-eta)*tand(theta)/h;
    end
end
% Calculate delta_ksi

% Initizalize vectors
F1 = zeros(1,NY);
F2 = zeros(1,NY);
F3 = zeros(1,NY);
F4 = zeros(1,NY);
F1_bar = zeros(1,NY);
F2_bar = zeros(1,NY);
F3_bar = zeros(1,NY);
F4_bar = zeros(1,NY);
G1 = zeros(1,NY);
G2 = zeros(1,NY);
G3 = zeros(1,NY);
G4 = zeros(1,NY);
df_F1_ks = zeros(1,NY);
df_F2_ks = zeros(1,NY);
df_F3_ks = zeros(1,NY);
df_F4_ks = zeros(1,NY);
u_mat = zeros(iter,NY);
v_mat = zeros(iter,NY);
T_mat = zeros(iter,NY);
p_mat = zeros(iter,NY);
rho_mat = zeros(iter,NY);
a_mat = zeros(iter,NY);
M_mat = zeros(iter,NY);

% Obtain initial conditions
F1(1) = rho1*u1;
F2(1) = rho1*u1^2 + p1;
F3(1) = rho1*u1*v1;
F4(1) = gamma/(gamma-1)*p1*u1 + rho1*u1*(u1^2+v1^2)/2;
G1(1) = rho1*v1;
G1(2) = rho1*u1*v1;
G1(3) = rho1*v1^2 + p1;
G1(4) = gamma/(gamma-1)*p1*v1 + rho1*u1*(u1^2+v1^2)/2;

% Main loop
for loop = 1:iter
    
    % Obtain geometry values
    h = calc_h(x,E,H,theta);
    ys = calc_ys(x,E,theta);
    df_eta_x = calc_df_eta_x(x,E,eta,h,theta);
    
    % Obtain forward step
    delta_ksi = 
    
    % First part
    for i = 2:NY-1
    
        % Predictor step
        df_F1_ks(i) = df_eta_x*(F1(i) - F1(i+1))/delta_eta + 1/h*(G1(i)-G1(i+1))/delta_eta;
        df_F2_ks(i) = df_eta_x*(F2(i) - F2(i+1))/delta_eta + 1/h*(G2(i)-G2(i+1))/delta_eta;
        df_F3_ks(i) = df_eta_x*(F3(i) - F3(i+1))/delta_eta + 1/h*(G3(i)-G3(i+1))/delta_eta;
        df_F4_ks(i) = df_eta_x*(F4(i) - F4(i+1))/delta_eta + 1/h*(G4(i)-G4(i+1))/delta_eta;
        
        % Obtain artificial viscosity terms
        SF1 = C*abs(p_mat(loop,i+1)-2*p_mat(loop,i)+p_mat(loop,i-1))/...
              (p_mat(loop,i+1)+2*p_mat(loop,i)-p_mat(loop,i-1)))*(F1(i+1)-2*F1(i)+F1(i-1));
        SF2 = C*abs(p_mat(loop,i+1)-2*p_mat(loop,i)+p_mat(loop,i-1))/...
              (p_mat(loop,i+1)+2*p_mat(loop,i)-p_mat(loop,i-1)))*(F2(i+1)-2*F2(i)+F2(i-1));
        SF3 = C*abs(p_mat(loop,i+1)-2*p_mat(loop,i)+p_mat(loop,i-1))/...
              (p_mat(loop,i+1)+2*p_mat(loop,i)-p_mat(loop,i-1)))*(F3(i+1)-2*F3(i)+F3(i-1));
        SF4 = C*abs(p_mat(loop,i+1)-2*p_mat(loop,i)+p_mat(loop,i-1))/...
              (p_mat(loop,i+1)+2*p_mat(loop,i)-p_mat(loop,i-1)))*(F4(i+1)-2*F4(i)+F4(i-1));

        % Calculate estimated values
        F1_bar(i) = F1(i) + df_F1_ks(i)*delta_ksi + SF1;
        F2_bar(i) = F1(i) + df_F2_ks(i)*delta_ksi + SF2;
        F3_bar(i) = F1(i) + df_F3_ks(i)*delta_ksi + SF3;
        F4_bar(i) = F1(i) + df_F4_ks(i)*delta_ksi + SF4;
    
    end    
    
    % Obtain estimated density
    A = F3.^2./(2*F1) - F4;
    B = gamma/(gamma-1)*F1.*F2;
    C = -(gamma+1)/(2*(gamma-1))*F1.^3;
    rho_bar = (-B+sqrt(B.^2-4*A.*C))./(2*A);
    
    % Obtain estimated pressure
    p_bar = F2_bar - F1_bar^2./rho_bar;
    
    % Obtain estimated G
    G1_bar = rho_bar.*F3_bar./F1_bar;
    G2_bar = F3_bar;
    G3_bar = rho_bar.*(F3_bar./F1_bar).^2 + F2_bar - F1_bar^2./rho_bar;
    G4_bar = gamma/(gamma-1)(F2_bar-F1_bar.^2./rho_bar)*F3_bar./F1_bar +...
             rho_bar/2.*F3_bar./F1_bar.*((F1_bar./rho_bar).^2+(F3_bar./F1_bar).^2);
    
    % Second part
    for i = 2:NY-1
        
        % Corrector step
        df_F1_ks_bar = df_eta_x*(F1_bar(i-1) - F1_bar(i))/delta_eta ...
                       + 1/h*(G1_bar(i-1)-G1_bar(i))/delta_eta;
        df_F2_ks_bar = df_eta_x*(F2_bar(i-1) - F2_bar(i))/delta_eta ...
                       + 1/h*(G2_bar(i-1)-G2_bar(i))/delta_eta;
        df_F3_ks_bar = df_eta_x*(F3_bar(i-1) - F3_bar(i))/delta_eta ...
                       + 1/h*(G3_bar(i-1)-G3_bar(i))/delta_eta;
        df_F4_ks_bar = df_eta_x*(F4_bar(i-1) - F4_bar(i))/delta_eta ...
                       + 1/h*(G4_bar(i-1)-G4_bar(i))/delta_eta;
    
        % Obtain artificial viscosity terms
        SF1_bar = 
        
        % Average derivatives
        df_F1_AV = (df_F1_ks(i) + df_F1_ks_bar(i))*delta_ksi/2;
        df_F2_AV = (df_F2_ks(i) + df_F2_ks_bar(i))*delta_ksi/2;
        df_F3_AV = (df_F3_ks(i) + df_F3_ks_bar(i))*delta_ksi/2;
        df_F4_AV = (df_F4_ks(i) + df_F4_ks_bar(i))*delta_ksi/2;
        
        
    end
    
end








