clc,clear

% Euler Implícito

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

S = alpha*delta_t/delta_x^2;

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

for loop = 2:iter+1 % Queremos os dados a partir da segunda linha da matriz sol
    % Solução numérica
    % Montar matriz A
    A = zeros(Im-1);
    % Primeira linha
    A(1,1:2) = [(1 + 2*S),-S];
    % Linhas entre a primeira e a última
    for i = 2:Im-2
        A(i,i-1:i+1) = [-S,(1 + 2*S),-S];
    end
    % Útima linha
    A(end,end-1:end) = [-S,(1 + 2*S)];
    % Vetor r
    r = sol(loop-1,2:Im);
    r(1) = r(1) + S*sol(loop,1);
    r(end) = r(end) + S*sol(loop,end);
    x_thomas = thomas_algorithm(A,r);
    sol(loop,2:end-1) = x_thomas;
    
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