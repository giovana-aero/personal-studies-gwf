clc,clear

% Euler Explícito

% Uma barra de comprimento L e inicialmente a uma temperatura de T_init.
% A extremidade direita é aquecida para T_Im e a extremidade esquerda é 
% mantida em T_0.

% Propriedades físicas
L = 1; % [m]
T_init = 0; % [°C]
T_Im = 100; % [°C]
T_0 = 0; % [°C]
alpha = 0.0834; % [m^2/s]

% Dados da malha e da solução
Im = 5; % (seis pontos no total, cinco elementos)
delta_x = L/Im;
s = 1/20; % (cumpre o requisito s <= 1/2)
%delta_t = s*delta_x^2/alpha;
delta_t = 0.1;
iter = 4/delta_t; % Número de iterações

% i -> ponto na malha
% n -> ponto no tempo

% Definir equação principal do método
function T_i__np1 = euler_ex(T_i__n,T_im1__n,T_ip1__n,delta_x,delta_t,alpha)
    T_i__np1 = T_i__n + delta_t*alpha/delta_x^2*(T_im1__n - 2*T_i__n + T_ip1__n); 
end

% Matriz pra guardar as soluções
sol = zeros(iter+1,Im+1); % Os seis pontos na barra, em cada instante de tempo

% Já se sabe a temperatura da barra inteira em t_0:
sol(1,:) = repmat(T_init,1,Im+1);

% Também ja se sabe as temperaturas das pontas:
sol(:,1) = repmat(T_0,iter+1,1);
sol(:,end) = repmat(T_Im,iter+1,1);

% Obter solução analítica
sol_an = sol;
x = linspace(0,L,Im+1);
t = 0;
%sol_an = [T_0,zeros(1,Im-1),T_Im];
%for i = 2:Im
%    sol_an(i) = heated_bar_analytic(x,L,T_Im,alpha,t,sum_max)

for loop = 2:iter+1 % Queremos os dados a partir da segunda linha da matriz sol
    % Solução numérica
    for node = 2:Im % Calcula-se apenas os pontos entre o primeiro e o último
        sol(loop,node) = euler_ex(sol(loop-1,node),sol(loop-1,node-1),sol(loop-1,node+1),delta_x,delta_t,alpha);
    end
    
    % Solução analítica
    for node = 2:Im
        sol_an(loop,node) = heated_bar_analytic(x(node),L,T_Im,alpha,t,100);
    end
    
    % Determinar o próximo instante no tempo
    % (importante pra solução analítica)
    t = t + delta_t;
    
    figure(1),clf
    plot([1:(Im+1)],sol(loop,:),'r','linewidth',2),grid on,hold on
    scatter([1:(Im+1)],sol(loop,:),'r')
    plot([1:(Im+1)],sol_an(loop,:),'k--','linewidth',2),grid on,hold on
    scatter([1:(Im+1)],sol_an(loop,:),'k')
    title(strcat('Iteração ',' ',num2str(loop-1)))
    xlabel('Malha'),ylabel('Temperatura [°C]')
    pause(0.05)
end

disp(sol(end,:))
disp(sol_an(end,:))
