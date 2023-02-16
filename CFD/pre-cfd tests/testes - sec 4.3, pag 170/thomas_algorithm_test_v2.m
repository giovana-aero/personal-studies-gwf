clc,clear

%https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm

% Entradas: matrizes A e r
A = [1,2,0,0,0,0;
     3,4,5,0,0,0;
     0,6,3,4,0,0;
     0,0,9,1,2,0;
     0,0,0,3,4,5;
     0,0,0,0,1,2]; % Tamanho N
r = [10;3;5;5;0;20];
%A = [1.1,-0.05,0,0;
%     -0.05,1.1,-0.05,0;
%     0,-0.05,1.1,-0.05;
%     0,0,-0.05,1.1];
%r = [0;0;0;5];
%r = [0.0004;0.0094;0.2074;9.5548];

disp(linsolve(A,r)')


N = length(r);
C = zeros(1,N);
D = zeros(1,N);

C(1) = A(1,2)/A(1,1);
D(1) = r(1)/A(1,1);
for i = 2:N
    if i < N
        C(i) = A(i,i+1)/(A(i,i) - A(i,i-1)*C(i-1));
    end
    
    D(i) = (r(i) - A(i,i-1)*D(i-1))/(A(i,i) - A(i,i-1)*C(i-1));
end

x = zeros(1,N);
x(N) = D(N);
for i = N-1:-1:1
    x(i) = D(i) - C(i)*x(i+1);
end
disp(x)