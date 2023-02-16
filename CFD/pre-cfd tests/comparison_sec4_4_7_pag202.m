clc,clear

% Dados do problema
L = 1;

% Condições de contorno
%val_down = 0;
%val_up = x^2;
%val_left = 0;
%val_right = y^2;

% Dados da malha
Im = 4;
Jm = 4;
delta_x = L/Im;
delta_y = L/Jm;
beta = delta_x/delta_y;

% Configurações da simulação
iter = 500;
ep = 1e-4;

% Solução analítica
an_sol = @(x,y) x^2*y^3;

% Montar vetores com coordenadas dos pontos
coo_x = linspace(0,L,Im+1);
coo_y = linspace(0,L,Jm+1);

% Inicializar primeira matriz de solução
sol1 = zeros(Jm+1,Im+1);

% Atribuir condições de contorno
sol1(end,:) = coo_x.^2;
sol1(:,end) = coo_y'.^2;

% Inicializar as demais matrizes de solução
sol2 = sol1;
sol3 = sol1;
sol4 = sol1;

% Gauss por ponto
disp("Gauss por ponto")
% Guardar a solução antes da próxima iteração pra comparar
sol_old = sol1;
for loop1 = 1:iter
    disp(["Iteração ",num2str(loop1)])
    for j = 2:Jm
        for i = 2:Im
            sol1(i,j) = (2*coo_y(j)*delta_x^2*(3*coo_x(i)^2 + coo_y(j)^2) - sol1(j,i+1) - sol1(j,i-1) - beta^2*sol1(j+1,i) - beta^2*sol1(j-1,i))*(-2 - 2*beta^2);
        end
    end
    
    % Checar convergência 
    sum_val = 0;
    for i = 2:Im
        for j = 2:Jm
            sum_val = sum_val + abs(sol1(j,i) - sol_old(j,i));
        end
    end
    
    if sum_val <= ep
        disp('Convergência')
        break
    end
    
    % Substituir pra fazer a devida comparação na iteração seguinte
    sol_old = sol1;
end

% Gauss por linha


% SOR por ponto


% SOR por linha
