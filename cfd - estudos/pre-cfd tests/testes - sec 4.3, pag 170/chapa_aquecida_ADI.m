clc,clear

% NOTA: o livro usao padrão (coluna,linha) ao invés de (linha,coluna).
% fazer as correções devidas depois

% Uma chapa aquecida quadrada de lados de comprimento L, difusividade térmica
% alpha, temperatura inicial T_init e temperaturas especificadas em cada lado

% Metodo ADI (alternating-direction implicit)

% Dados do problema
L = 1; % [m]
alpha = 0.0155; % [m^2/s]
T_init = 25; % [°C]
T_up = 100; % [°C]
T_down = 20; % [°C]
T_right = 30; % [°C]
T_left = 30; % [°C]

% Dados da malha e da solução
Im = 5; % horizontal (número total de pontos = Im + 1)
Jm = 5; % vertical (número total de pontos = Jm + 1)
delta_x = L/Im;
delta_y = L/Jm;
delta_t = 1;
sx = alpha*delta_t/delta_x^2; S = sx;
%sy = alpha*delta_t/delta_y^2;
iter = 20/delta_t;

% Matriz da solução:
sol = zeros(Jm+1,Im+1) + T_init; % Já deixando com a temperatura inicial

% Também já se sabe quais são as condições de contorno:
sol(1,:) = repmat(T_up,1,Im+1);
sol(end,:) = repmat(T_down,1,Im+1);
sol(:,1) = repmat(T_left,Jm+1,1);
sol(:,end) = repmat(T_right,Jm+1,1);

% Solução
for loop = 1:iter
    % Montar matriz A
    A = zeros(Im-1);
    % Primeira linha
    A(1,1:2) = [2*(1 + S),-S];
    % Linhas entre a primeira e a última
    for i = 2:Im-2
        A(i,i-1:i+1) = [-S,2*(1 + S),-S];
    end
    % Útima linha
    A(end,end-1:end) = [-S,2*(1 + S)];
    
    % Montar vetor r
    temp = flip(sol); % Fazer uma matriz temporária pra facilitar
    r = S*temp(3,2:Jm) + 2*(1-S)*temp(2,2:Jm) + S*temp(1,2:Jm);
    r(1) = r(1) + S*temp(2,1);
    r(end) = r(end) + S*temp(2,end);
    
    % Solucionar 
    x = thomas_algorithm(A,r)
    
    
    
    
end