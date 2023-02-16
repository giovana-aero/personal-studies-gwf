clc,clear

% Método de Gauss-Seidel por Ponto (PGS)
% Condição de fronteira tipo Neumann na fronteira direita:
% Fluxo de calor nulo

% Dados do problema
T_init = 0; % [°C]
T_up = 0; % [°C]
T_down = 50; % [°C]
%T_right = 100; % [°C]
T_left = 75; % [°C]

% Dados da geometria e da malha
L = 1; % [m] (chapa quadrada)
Im = 4;
Jm = 4;
delta_x = L/(Im+1);
delta_y = L/(Jm+1);
beta = delta_x/delta_y;
x_dots = linspace(0,L,Im+1);
y_dots = linspace(0,L,Jm+1);

% Iterações e critério de parada
iter = 100;
ep = 0.0001;

% Criar matriz de solução
sol = zeros(Jm+1,Im+2) + T_init; % Adicionar mais uma coluna pra representar os pontos fantasmas à direita

% Atribuir as condições de contorno
sol(1,:) = T_up;
sol(end,:) = T_down;
%sol(:,end) = T_right;
sol(:,1) = T_left;

% Virar de cabeça pra baixo pra facilitar o processo, considerando
% o padrão adotado pelo livro
sol = flip(sol);
% Guardar a solução antes da próxima iteração pra comparar
sol_old = sol;

for loop = 1:iter
    disp(["Iteração ",num2str(loop)])
    
    % Cumprir a condição da parede isolada
%    sol(:,end) = sol(:,end-2);
    
    % Avaliar as equações nos pontos internos da matriz
    for j = 2:Jm 
        for i = 2:Im+1
            if i < Im+1
                sol(j,i) = 1/(2*(1+beta^2))*(sol(j,i+1) + sol(j,i-1) + beta^2*sol(j+1,i) + beta^2*sol(j-1,i));
            else
                sol(j,i) = 1/(2*(1+beta^2))*(2*sol(j,i-1) + beta^2*sol(j+1,i) + beta^2*sol(j-1,i));
            end
            if i == Im
                sol(j,end) = sol(j,i);
            end
        end
    end
    
    % Fazer um gráfico (desconsiderar os pontos fantasmas)
    figure(1),clf
    contourf(x_dots,y_dots,flip(sol(:,1:end-1))),grid on
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
disp(flip(sol(:,1:end-1)))