clc,clear

% Método de Gauss-Seidel por Ponto (PGS)

% Dados do problema
T_init = 0; % [°C]
T_up = 0; % [°C]
T_down = 50; % [°C]
T_right = 100; % [°C]
T_left = 75; % [°C]

% Dados da geometria e da malha
L = 1; % [m] (chapa quadrada)
Im = 20;
Jm = 20;
delta_x = L/(Im+1);
delta_y = L/(Jm+1);
beta = delta_x/delta_y;
x_dots = linspace(0,L,Im+1);
y_dots = linspace(0,L,Jm+1);

% Iterações e critério de parada
iter = 1000;
ep = 0.0001;

% Criar matriz de solução
sol = zeros(Jm+1,Im+1) + T_init;

% Atribuir as condições de contorno
sol(1,:) = T_up;
sol(end,:) = T_down;
sol(:,end) = T_right;
sol(:,1) = T_left;

% Virar de cabeça pra baixo pra facilitar o processo, considerando
% o padrão adotado pelo livro
sol = flip(sol);
% Guardar a solução antes da próxima iteração pra comparar
sol_old = sol;

for loop = 1:iter
    disp(["Iteração ",num2str(loop)])
    % Avaliar as equações nos pontos internos da matriz
    for j = 2:Jm 
        for i = 2:Im
            sol(j,i) = 1/(2*(1+beta^2))*(sol(j,i+1) + sol(j,i-1) + beta^2*sol(j+1,i) + beta^2*sol(j-1,i));
        end
    end
    
    % Fazer um gráfico
    figure(1),clf
    contourf(x_dots,y_dots,flip(sol)),grid on
    xlabel('x'),ylabel('y')
    title(['Iteração ', num2str(loop)])
    c = colorbar;
    ylabel(c,"Temperatura [°C]")
    pause(0.001)
    
    % Checar convergência 
    sum_val = 0;
    for i = 2:Im
        for j = 2:Jm
            sum_val = sum_val + abs(sol(j,i) - sol_old(j,i));
        end
    end
    
    if sum_val <= ep
        disp('Convergência')
        break
    end
    
    % Substituir pra fazer a devida comparação na iteração seguinte
    sol_old = sol;
    
end

disp('')
disp('Solução:')
disp(flip(sol))