clc,clear

% Entradas: matrizes A e r
%A = [1,2,0,0,0,0;
%     3,4,5,0,0,0;
%     0,6,3,4,0,0;
%     0,0,9,1,2,0;
%     0,0,0,3,4,5;
%     0,0,0,0,1,2]; % Tamanho N
%r = [10;3;5;5;0;20];
A = [1.1,-0.05,0,0;
     -0.05,1.1,-0.05,0;
     0,-0.05,1.1,-0.05;
     0,0,-0.05,1.1];
r = [0;0;0;5];
%r = [0.0004;0.0094;0.2074;9.5548];

disp(linsolve(A,r)')

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
disp(x)