function x = thomas_algorithm(A,r)

% Entradas: matrizes A e r


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