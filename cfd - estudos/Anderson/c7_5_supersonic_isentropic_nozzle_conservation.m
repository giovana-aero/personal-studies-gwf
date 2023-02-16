clc,clear

% Solution data
iter = 1500; % Number of iterations (main loop)
C = 0.5; % Courant number
gamma = 1.4;
    
% Mesh data
N = 31; % Number of points
x_vec = linspace(0,3,N); % Distribution of points
delta_x = 3/(N-1);

% Nozze shape vector
A_vec = 1 + 2.2*(x_vec - 1.5).^2;
[~,TH_pos] = min(A_vec); % Throat position (or the closest point)

% Initialize solution vectors (all dimensionless)
rho_vec = zeros(1,N);
T_vec = zeros(1,N);
for i = 1:length(x_vec)
    if x_vec(i) <= 0.5
        rho_vec(i) = 1;
        T_vec(i) = 1;
    elseif x_vec(i) > 0.5 && x_vec(i) <= 1.5
        rho_vec(i) = 1 - 0.366*(x_vec(i)-0.5);
        T_vec(i) = 1 - 0.167*(x_vec(i)-0.5);
    else
        rho_vec(i) = 0.634 - 0.3879*(x_vec(i)-1.5);
        T_vec(i) = 0.833 - 0.3507*(x_vec(i)-1.5);
    end
end
a_vec = sqrt(T_vec);
V_vec = 0.59./(rho_vec.*A_vec);
mach_vec = V_vec./a_vec;
U1 = rho_vec.*A_vec;
U2 = rho_vec.*A_vec.*V_vec;
U3 = rho_vec.*(T_vec/(gamma-1)+gamma/2*V_vec.^2).*A_vec;

% Calculate delta_t
function ans = calc_delta_t(C,delta_x,a_vec,V_vec)
    % Take the minimum value of all possible time steps
    vec = C*delta_x./(a_vec + V_vec);
    ans = min(vec);
end

% Vectors for historics (throat)
%his_U1 = zeros(1,iter);
%his_U2 = zeros(1,iter);
%his_U3 = zeros(1,iter);
his_V = zeros(1,iter);
his_mach = zeros(1,iter);
his_rho = zeros(1,iter);
his_T = zeros(1,iter);
his_p = zeros(1,iter);
his_res = zeros(3,iter);

% Main loop
for loop = 1:iter
    disp(['Iteration ',num2str(loop)])
    
    % Obtain time step
    delta_t = calc_delta_t(C,delta_x,a_vec,V_vec);
    
    % Calculate predicted values at point 1    
    % Calculate relevant terms
    F1i = U2(1);
    F1ip = U2(2);
    F2i = U2(1)^2/U1(1) + (gamma-1)/gamma*(U3(1)-gamma/2*U2(1)^2/U1(1));
    F2ip = U2(2)^2/U1(2) + (gamma-1)/gamma*(U3(2)-gamma/2*U2(2)^2/U1(2));
    F3i = gamma*U2(1)*U3(1)/U1(1) - gamma*(gamma-1)/2*U2(1)^3/U1(1)^2;
    F3ip = gamma*U2(2)*U3(2)/U1(2) - gamma*(gamma-1)/2*U2(2)^3/U1(2)^2;
    J2F = (gamma-1)/gamma*(U3(1)-gamma/2*U2(1)^2/U1(1))*(log(A_vec(2))-log(A_vec(1)))/delta_x;
    
    % Forward-differences
    diff_U1_F = -(F1ip - F1i)/delta_x;
    diff_U2_F = -(F2ip - F2i)/delta_x + J2F;
    diff_U3_F = -(F3ip - F3i)/delta_x;
    
    % Obtain predicted values
    U1_bar = U1(1) + diff_U1_F*delta_t;
    U2_bar = U2(1) + diff_U2_F*delta_t;
    U3_bar = U3(1) + diff_U3_F*delta_t;

    % Keep variable values
    % Calculate relevant terms
    F1_bar_old = U2_bar;
    F2_bar_old = U2_bar^2/U1_bar + (gamma-1)/gamma*(U3_bar-gamma/2*U2_bar^2/U1_bar);
    F3_bar_old = gamma*U2_bar*U3_bar/U1_bar - gamma*(gamma-1)/2*U2_bar^3/U1_bar^2;
    
    % Keep variable values
%    rho_vec_old = rho_vec;
%    T_vec_old = T_vec;
%    V_vec_old = V_vec;
%    a_vec_old = a_vec;
    
    % Calculate properties of all internal points
    for i = 2:N-1
        % Calculate relevant terms
        F1i = U2(i);
        F1ip = U2(i+1);
        F2i = U2(i)^2/U1(i) + (gamma-1)/gamma*(U3(i)-gamma/2*U2(i)^2/U1(i));
        F2ip = U2(i+1)^2/U1(i+1) + (gamma-1)/gamma*(U3(i+1)-gamma/2*U2(i+1)^2/U1(i+1));
        F3i = gamma*U2(i)*U3(i)/U1(i) - gamma*(gamma-1)/2*U2(i)^3/U1(i)^2;
        F3ip = gamma*U2(i+1)*U3(i+1)/U1(i+1) - gamma*(gamma-1)/2*U2(i+1)^3/U1(i+1)^2;
        J2F = (gamma-1)/gamma*(U3(i)-gamma/2*U2(i)^2/U1(i))*(log(A_vec(i+1))-log(A_vec(i)))/delta_x;
        
        % Forward-differences
        diff_U1_F = -(F1ip - F1i)/delta_x;
        diff_U2_F = -(F2ip - F2i)/delta_x + J2F;
        diff_U3_F = -(F3ip - F3i)/delta_x;
        
        % Obtain predicted values
        U1_bar = U1(i) + diff_U1_F*delta_t;
        U2_bar = U2(i) + diff_U2_F*delta_t;
        U3_bar = U3(i) + diff_U3_F*delta_t;
        
        % Calculate relevant terms
        F1_bar = U2_bar;
        F2_bar = U2_bar^2/U1_bar + (gamma-1)/gamma*(U3_bar-gamma/2*U2_bar^2/U1_bar);
        F3_bar = gamma*U2_bar*U3_bar/U1_bar - gamma*(gamma-1)/2*U2_bar^3/U1_bar^2;
        J2R = (gamma-1)/gamma*(U3_bar-gamma/2*U2_bar^2/U1_bar)*(log(A_vec(i))-log(A_vec(i-1)))/delta_x;
        
        % Rearward-differences
        diff_U1_R = -(F1_bar - F1_bar_old)/delta_x;
        diff_U2_R = -(F2_bar - F2_bar_old)/delta_x + J2R;
        diff_U3_R = -(F3_bar - F3_bar_old)/delta_x;
           
        % Average time derivatives
        diff_U1_AV = (diff_U1_F + diff_U1_R)/2;
        diff_U2_AV = (diff_U2_F + diff_U2_R)/2;
        diff_U3_AV = (diff_U3_F + diff_U3_R)/2;
        
        % Corrected values
        U1(i) = U1(i) + diff_U1_AV*delta_t;
        U2(i) = U2(i) + diff_U2_AV*delta_t;
        U3(i) = U3(i) + diff_U3_AV*delta_t;
        
        % Keep predicted values to use at the next point
        F1_bar_old = F1_bar;
        F2_bar_old = F2_bar;
        F3_bar_old = F3_bar;
%        U1_bar_old = U1_bar;
%        U2_bar_old = U2_bar;
%        U3_bar_old = U3_bar;
        
        % Keep residuals at the throat
        if i == TH_pos
            his_res(:,loop) = abs([diff_U1_AV;diff_U2_AV;diff_U3_AV]);
        end
    end
    
    % Recalculate all flow properties
    rho_vec = U1./A_vec;
    V_vec = U2./U1;
    T_vec = (gamma-1)*(U3./U1-gamma/2*V_vec.^2);
    p_vec = rho_vec.*T_vec;
    a_vec = sqrt(T_vec);
    mach_vec = V_vec./a_vec;
    
    % Apply boundary conditions
    U2(1) = 2*U2(2) - U2(3);
    U3(1) = U1(1)*(T_vec(1)/(gamma-1)+gamma/2*V_vec(1)^2);
    U1(N) = 2*U1(N-1) - U1(N-2);
    U2(N) = 2*U2(N-1) - U2(N-2);
    U3(N) = 2*U3(N-1) - U3(N-2);
    
    % Store values at the throat
    his_rho(loop) = rho_vec(TH_pos);
    his_T(loop) = T_vec(TH_pos);
    his_p(loop) = p_vec(TH_pos);
    his_V(loop) = V_vec(TH_pos);
    his_mach(loop) = mach_vec(TH_pos);
    
end




% Plot a graph
figure(1),clf
plot(x_vec,A_vec/A_vec(1)),hold on,grid on
plot(x_vec,V_vec)
plot(x_vec,mach_vec)
plot(x_vec,rho_vec)
plot(x_vec,T_vec)
plot(x_vec,p_vec)
%title(['Iteration ',num2str(loop)])
legend('Nozzle contour','Velocity','Mach','Density','Temperature','Pressure')
title('Properties along the complete nozzle')

% Plot properties at the throat
figure(2),clf
plot([1:iter],his_V),hold on,grid on
plot([1:iter],his_mach)
plot([1:iter],his_rho)
plot([1:iter],his_T)
plot([1:iter],his_p)
legend('Velocity','Mach','Density','Temperature','Pressure')
title('Properties at the throat')

% Plot residuals
figure(3),clf
plot([1:iter],his_res(1,:)),hold on,grid on
plot([1:iter],his_res(2,:))
plot([1:iter],his_res(3,:))
legend('Velocity','Density','Temperature')
title('Residuals (absolute values)')