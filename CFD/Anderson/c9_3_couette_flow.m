clc,clear

% Mesh data
N = 21;
delta_y = 1/(N-1);

% Problem data
Re = 5000;
E = 1; % Parameter for the calculation of delta_t
iter = 10;
conv = eps; % Convergence criterion

% Initialize solution
u_vec = zeros(1,N);
u_vec(N) = 1;


% Functions
function [A,B] = calc_A_B(delta_t,delta_y,Re)
    A = -delta_t/(2*delta_y^2*Re);
    B = 1 + delta_t/(delta_y^2*Re);
end
function K = calc_K(delta_t,delta_y,Re,uj,ujp1,ujm1)
    K = (1-delta_t/(Re*delta_y^2))*uj + delta_t/(2*delta_y^2*Re)*(ujp1+ujm1);
end
function delta_t = calc_delta_t(Re,delta_y,E)
    delta_t = E*Re*delta_y^2;
end
function x = thomas_algorithm(A,r)
    % Determinar os valores de alpha, beta e gamma
    N = length(r);
    alpha = zeros(1,N);
    beta = zeros(1,N);
    gamma = zeros(1,N);
    % Passo 1 (i=1)
    beta(1) = A(1,1);
    gamma(1) = A(1,2)/beta(1);
    % Passos 2 e 3 (2 <= i <= N)
    for i = 2:N
        alpha(i) = A(i,i-1);
        beta(i) = A(i,i) - A(i,i-1)*gamma(i-1);
        if i < N
            gamma(i) = A(i,i+1)/beta(i);
        end
    end
    % Determinar os valores de si
    s = zeros(1,N);
    s(1) = r(1)/beta(1);
    for i = 2:N
        s(i) = (r(i) - alpha(i)*s(i-1))/beta(i);
    end
    % Obter vetor solução x
    x = zeros(1,N);
    x(N) = s(N);
    for i = N-1:-1:1
        x(i) = s(i) - gamma(i)*x(i+1);
    end
end

% Obtain step time
delta_t = calc_delta_t(Re,delta_y,E);

% Main loop
for i = 1:iter
    disp(['Iteration ',num2str(i)])

    % Initialize solution matrices
    mat1 = zeros(N-2);
    mat2 = zeros(N-2,1);

    % Obtain A & B
    [A,B] = calc_A_B(delta_t,delta_y,Re);

    for j = 2:N-1
        % Obtain K
        K = calc_K(delta_t,delta_y,Re,u_vec(j),u_vec(j+1),u_vec(j-1));

        % Insert into solution matrices
        if j == 2
            mat1(1,1) = B;
            mat1(1,2) = A;
            mat2(j-1) = K;
        elseif j == N-1
            mat1(end,end) = B;
            mat1(end,end-1) = A;
            mat2(j-1) = K - A*1;
        else
            mat1(j-1,j-1) = B;
            mat1(j-1,j) = A;
            mat1(j-1,j-2) = A;
            mat2(j-1) = K;
        end
    end

    % Keep old values to compare
    u_old = u_vec;

    % Solve system using Thomas' algorithm
    u_vec = thomas_algorithm(mat1,mat2);
    u_vec = [0,u_vec,1];

    % Check for convergence
    vec = u_vec-u_old;
    if max(abs(vec)) <= eps
        disp('Convergence')
        break
    end
end

% Plot result
%figure(1),clf
%plot([1:N],u_vec,'k'),grid on

figure(2),clf
plot(u_vec,linspace(0,1,N),'k'),grid on
xlabel('u/ue'),ylabel('y/D')
