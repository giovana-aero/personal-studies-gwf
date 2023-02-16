clc,clear

% Solution data
iter = 3000; % Number of iterations (main loop)
C = 0.4; % Courant number
gamma = 1.4;
p_N = 0.93; % Exit pressure ratio
    
% Mesh data
N = 31; % Number of points
x_vec = linspace(0,3,N); % Distribution of points
delta_x = 3/(N-1);

% Initialize solution vectors (all dimensionless)
rho_vec = 1 - 0.023*x_vec;
T_vec = 1 - 0.009333*x_vec;
V_vec = 0.05 + 0.11*x_vec;
%rho_vec = ones(1,N);
%T_vec = ones(1,N);
%V_vec = ones(1,N);
a_vec = sqrt(T_vec);
p_vec = ones(1,N); p(N) = p_N;
mach_vec = zeros(1,N);

% Nozze shape vector
A_vec = zeros(1,N);
for i = 1:length(x_vec)
    if x_vec(i) <= 1.5
        A_vec(i) = 1 + 2.2*(x_vec(i)-1.5)^2;
    else
        A_vec(i) = 1 + 0.2223*(x_vec(i)-1.5)^2;
    end
end
%A_vec = 1 + 2.2*(x_vec - 1.5).^2;
[~,TH_pos] = min(A_vec); % Throat position (or the closest point)

% Calculate delta_t
function ans = calc_delta_t(C,delta_x,a_vec,V_vec)
    % Take the minimum value of all possible time steps
    vec = C*delta_x./(a_vec + V_vec);
    ans = min(vec);
end

% Vectors for historics (throat)
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
    % Forward-differences
    diff_rho_F = -rho_vec(1)*(V_vec(2) - V_vec(1))/delta_x - ...
                 rho_vec(1)*V_vec(1)*(log(A_vec(2))-log(A_vec(1)))/delta_x - ...
                 V_vec(1)*(rho_vec(2)-rho_vec(1))/delta_x;
    diff_V_F = -V_vec(1)*(V_vec(2)-V_vec(1))/delta_x - 1/gamma*((T_vec(2)-T_vec(1))/delta_x + ...
                T_vec(1)/rho_vec(1)*(rho_vec(2)-rho_vec(1))/delta_x);
    diff_T_F = -V_vec(1)*(T_vec(2)-T_vec(1))/delta_x - (gamma-1)*T_vec(1)*((V_vec(2)-V_vec(1))/delta_x + ...
                V_vec(1)*(log(A_vec(2))-log(A_vec(1)))/delta_x);
    % Obtain predicted values
    rho_bar_old = rho_vec(1) + diff_rho_F*delta_t;
    V_bar_old = V_vec(1) + diff_V_F*delta_t;
    T_bar_old = T_vec(1) + diff_T_F*delta_t;
    
    % Keep variable values
%    rho_vec_old = rho_vec;
%    T_vec_old = T_vec;
%    V_vec_old = V_vec;
%    a_vec_old = a_vec;
    
    % Calculate properties of all internal points
    for i = 2:N-1
        % Forward-differences
        diff_rho_F = -rho_vec(i)*(V_vec(i+1) - V_vec(i))/delta_x - ...
                     rho_vec(i)*V_vec(i)*(log(A_vec(i+1))-log(A_vec(i)))/delta_x - ...
                     V_vec(i)*(rho_vec(i+1)-rho_vec(i))/delta_x;
        diff_V_F = -V_vec(i)*(V_vec(i+1)-V_vec(i))/delta_x - 1/gamma*((T_vec(i+1)-T_vec(i))/delta_x + ...
                    T_vec(i)/rho_vec(i)*(rho_vec(i+1)-rho_vec(i))/delta_x);
        diff_T_F = -V_vec(i)*(T_vec(i+1)-T_vec(i))/delta_x - (gamma-1)*T_vec(i)*((V_vec(i+1)-V_vec(i))/delta_x + ...
                    V_vec(i)*(log(A_vec(i+1))-log(A_vec(i)))/delta_x);
        % Obtain predicted values
        rho_bar = rho_vec(i) + diff_rho_F*delta_t;
        V_bar = V_vec(i) + diff_V_F*delta_t;
        T_bar = T_vec(i) + diff_T_F*delta_t;
        % Rearward-differences
        diff_rho_R = -rho_bar*(V_bar-V_bar_old)/delta_x - ...
                     rho_bar*V_bar*(log(A_vec(i))-log(A_vec(i-1)))/delta_x - ...
                     V_bar*(rho_bar-rho_bar_old)/delta_x;
        diff_V_R = -V_bar*(V_bar-V_bar_old)/delta_x - 1/gamma*((T_bar-T_bar_old)/delta_x + ...
                    T_bar/rho_bar*(rho_bar-rho_bar_old)/delta_x);
        diff_T_R = -V_bar*(T_bar-T_bar_old)/delta_x - (gamma-1)*T_vec(i)*((V_bar-V_bar_old)/delta_x + ...
                    V_bar*(log(A_vec(i))-log(A_vec(i-1)))/delta_x);        
        % Average time derivatives
        diff_rho_AV = (diff_rho_F + diff_rho_R)/2;
        diff_V_AV = (diff_V_F + diff_V_R)/2;
        diff_T_AV = (diff_T_F + diff_T_R)/2;
        % Corrected values
        rho_vec(i) = rho_vec(i) + diff_rho_AV*delta_t;
        V_vec(i) = V_vec(i) + diff_V_AV*delta_t;
        T_vec(i) = T_vec(i) + diff_T_AV*delta_t;
        p_vec(i) = rho_vec(i)*T_vec(i);
        % Keep predicted values to use at the next point
        rho_bar_old = rho_bar;
        V_bar_old = V_bar;
        T_bar_old = T_bar;
        % Keep residuals at the throat
        if i == TH_pos
            his_res(:,loop) = abs([diff_V_AV;diff_rho_AV;diff_T_AV]);
        end
    end
    
    % Calculate pressure at the end
%    p_vec(N) = rho_vec(N)*T_vec(N);
    
    % Recalculate a_vec
    a_vec = sqrt(T_vec);
    
    % Obtain Mach numbers
    mach_vec = V_vec./a_vec;
    
    % Apply boundary conditions
    V_vec(1) = 2*V_vec(2) - V_vec(3);
    V_vec(N) = 2*V_vec(N-1) - V_vec(N-2);
    T_vec(N) = 2*T_vec(N-1) - T_vec(N-2);
    rho_vec(N) = p_vec(N)/T_vec(N);
    
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