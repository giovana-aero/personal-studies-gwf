clc,clear

% Seção 6.2.1 - Escoamento em cavidade com MAC

% Assumindo que o fluido é água a 20°C

% Dados do problema
H = 1; % [m] altura
L = 1; % [m] largura
% Água a 20 °C
visc = 0.000001; % [m^2/s] (viscosidade cinemática) https://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html
rho = 1e3; % [kg/m^3]
% Ar no nível do mar
%visc = 14.6e-6; % retirado do mesmo site referenciado acima
%rho = 1.225; % Retirado da tabela A do fundamentos de engenharia aeronáutica
Re = 100;
u0 = Re*visc/H; % Velocidade da tampa

% Dados da malha e da solução
NI = 10; % Número de pontos horizontais
NJ = 10; % Número de pontos verticais
iter1 = 50; % Iterações pros cálculos de pressão
iter2 = 10; % Iterações da simulação completa
delta_x = L/(NI-2);
delta_y = H/(NJ-2);
w = 1.9; % Fator de relaxação do PSOR (checar página 197)
beta = delta_x/delta_y;

% Inicializar matrizes com resultados
u_mat = zeros(NJ,NI);
v_mat = zeros(NJ,NI);
p_mat = zeros(NJ,NI);
F_mat = zeros(NJ,NI);
G_mat = zeros(NJ,NI);

% Funções necessárias pras interpolações
function ans = sinal(x)
    if x >= 0
        ans = 1;
    else
        ans = -1;
    end
end
function ans = calc_pe(ui,delta_x,visc)
    ans = abs(ui)*delta_x/visc;
end
function ans = upwind_1st(ue,u1,u2)
    ans = (1 + sinal(ue))/2*u1 + (1 - sinal(ue))/2*u2;
end
function ans = central_diff(u1,u2)
    ans = (u1 + u2)/2;
end
function ans = hybrid_itp(pe,ue,u1,u2)
    if pe < 1.9
        ans = central_diff(u1,u2);
    elseif 1.9 <= pe && pe < 2
        FP = (pe - 1.9)/0.1;
        ans = (1 - FP)*central_diff(u1,u2) + FP*upwind_1st(ue,u1,u2);
    else
        ans = upwind_1st(ue,u1,u2);
    end
end
% Funções pra calcular CONV e VISC de F e G, respectivamente
function ans = CONV_F(i,j,u_mat,v_mat,delta_x,delta_y,visc)

    % Termo del(u*u)/(del x)
    % Velocidades
    u_bar_ip3o2_j = (u_mat(j,i+2) + u_mat(j,i+1))/2;
    u_bar_ip1o2_j = (u_mat(j,i) + u_mat(j,i+1))/2;
    % Interpolação de u_itp_ip3o2_j
    ue = (u_mat(j,i+2) + u_mat(j,i+1))/2;
    pe = calc_pe(ue,delta_x,visc);
    u_itp_ip3o2_j= hybrid_itp(pe,ue,u_mat(j,i+2),u_mat(j,i+1));
    % Interpolação de u_itp_ip1o2_j
    ue = (u_mat(j,i+1) + u_mat(j,i))/2;
    pe = calc_pe(ue,delta_x,visc);
    u_itp_ip1o2_j= hybrid_itp(pe,ue,u_mat(j,i+1),u_mat(j,i));

    % Termo del(u*v)/(del y)
    % Velocidades
    v_bar_ip1o2_jp1 = (v_mat(j+1,i+1) + v_mat(j+1,i))/2;
    v_bar_ip1o2_j = (v_mat(j,i+1) + v_mat(j,i))/2;
    % Interpolação de u_itp_ip1_jp1o2
    ve = (v_mat(j+1,i+1) + v_mat(j+1,i))/2;
    pe = calc_pe(ve,delta_y,visc);
    u_itp_ip1_jp1o2 = hybrid_itp(pe,ve,u_mat(j,i+1),u_mat(j+1,i+1));
    % Interpolação de u_itp_ip1_jm1o2
    ve = (v_mat(j,i+1) + v_mat(j,i))/2;
    pe = calc_pe(ve,delta_y,visc);
    u_itp_ip1_jm1o2 = hybrid_itp(pe,ve,u_mat(j,i+1),u_mat(j-1,i+1));

    % Calcular
    ans = (u_bar_ip3o2_j*u_itp_ip3o2_j - u_bar_ip1o2_j*u_itp_ip1o2_j)/delta_x + ...
           (v_bar_ip1o2_jp1*u_itp_ip1_jp1o2 - v_bar_ip1o2_j*u_itp_ip1_jm1o2)/delta_y;

end
function ans = VISC_F(i,j,u_mat,delta_x,delta_y,visc)
    % Tomar valores de u
    u_i_j = u_mat(j,i);
    u_ip1_j = u_mat(j,i+1);
    u_ip2_j = u_mat(j,i+2);
    u_ip1_jp1 = u_mat(j+1,i+1);
    u_ip1_jm1 = u_mat(j-1,i+1);
    % Calcular
    ans = visc*(u_i_j - 2*u_ip1_j + u_ip2_j)/delta_x^2 + visc*(u_ip1_jp1 - 2*u_ip1_j + u_ip1_jm1)/delta_y^2;
end
function ans = CONV_G(i,j,u_mat,v_mat,delta_x,delta_y,visc)

    % Termo del(vv)/(del y)
    % Velocidades
    v_bar_i_jp3o2 = (v_mat(j+1,i) + v_mat(j+2,i))/2;
    v_bar_i_jp1o2 = (v_mat(j,i) + v_mat(j+1,i))/2;
    % Interpolação de v_itp_i_jp3o2
    ve = (v_mat(j+1,i) + v_mat(j+2,i))/2;
    pe = calc_pe(ve,delta_y,visc);
    v_itp_i_jp3o2 = hybrid_itp(pe,ve,v_mat(j+1,i),v_mat(j+2,i));
    % Interpolação de v_itp_i_jp1o2
    ve = (v_mat(j,i) + v_mat(j+1,i))/2;
    pe = calc_pe(ve,delta_y,visc);
    v_itp_i_jp1o2 = hybrid_itp(pe,ve,v_mat(j,i),v_mat(j+1,i));

    % Termo del(uv)/(del x)
    % Velocidades
    u_bar_ip1_jp1o2 = (u_mat(j,i+1) + u_mat(j+1,i+1))/2;
    u_bar_i_jp1o2 = (u_mat(j,i) + u_mat(j+1,i))/2;
    % Interpolação de v_itp_ip1o2_jp1
    ue = (u_mat(j,i+1) + u_mat(j+1,i+1))/2;
    pe = calc_pe(ue,delta_x,visc);
    v_itp_ip1o2_jp1 = hybrid_itp(pe,ue,v_mat(j+1,i),v_mat(j+1,i+1));
    % Interpolação de v_itp_im1o2_jp1
    ue = (u_mat(j,i) + u_mat(j+1,i))/2;
    pe = calc_pe(ue,delta_x,visc);
    v_itp_im1o2_jp1 = hybrid_itp(pe,ue,v_mat(j+1,i-1),v_mat(j+1,i));

    % Calcular
    ans = (v_bar_i_jp3o2*v_itp_i_jp3o2 - v_bar_i_jp1o2*v_itp_i_jp1o2)/delta_y + ...
           (u_bar_ip1_jp1o2*v_itp_ip1o2_jp1 - u_bar_i_jp1o2*v_itp_im1o2_jp1)/delta_x;

end
function ans = VISC_G(i,j,v_mat,delta_x,delta_y,visc)
    % Tomar valores de v
    v_i_j = v_mat(j,i);
    v_i_jp1 = v_mat(j+1,i);
    v_i_jp2 = v_mat(j+2,i);
    v_im1_jp1 = v_mat(j+1,i-1);
    v_ip1_jp1 = v_mat(j+1,i+1);
    % Calcular
    ans = visc*(v_i_j - 2*v_i_jp1 + v_i_jp2)/delta_y^2 + visc*(v_im1_jp1 - 2*v_i_jp1 + v_ip1_jp1)/delta_x^2;
end

% Função pra calcular o devido delta_t
tau = 0.4; % Fator de segurança
function ans = calc_delta_t(u_mat,v_mat,delta_x,delta_y,visc,tau)
    ans = tau*(abs(max(max(u_mat)))/delta_x + abs(max(max(v_mat)))/delta_y + 2*visc*(delta_x^-2 + delta_y^-2))^-1;
end

% Dados pros gráficos
[X,Y] = meshgrid(linspace(-L/2,L/2,NI-2),linspace(-H/2,H/2,NJ-2));

% Loop principal
for loop2 = 1:iter2
    disp(['Solucionador - Iteração ',num2str(loop2)])

    % Aplicar condições de fronteira das velocidades
    for i = 3:NI-1
        % Parede inferior
        u_mat(1,i) = -u_mat(2,i); % (eq 6.88, pag 343)
        % Parede superior
        u_mat(NJ,i) = 2*u0 - u_mat(NJ-1,i); % (eq 6.90, pag 343)
    end
    for j = 3:NJ-1
        % Paredes laterais
        v_mat(j,1) = -v_mat(j,2); % (eq 6.89, pag 343)
        v_mat(j,NI) = -v_mat(j,NI-1);
    end

    % Obter delta_t
    delta_t = calc_delta_t(u_mat,v_mat,delta_x,delta_y,visc,tau);

    % Calcular F
    for i = 2:NI-2
        for j = 2:NJ-1
            F_mat(j,i+1) = u_mat(j,i+1) + delta_t*(-CONV_F(i,j,u_mat,v_mat,delta_x,delta_y,visc) + VISC_F(i,j,u_mat,delta_x,delta_y,visc)); % (eq 6.56, pag 329)
        end
    end

    % Calcular G
    for i = 2:NI-1
        for j = 2:NJ-2
            G_mat(j+1,i) = v_mat(j+1,i) + delta_t*(-CONV_G(i,j,u_mat,v_mat,delta_x,delta_y,visc) + VISC_G(i,j,v_mat,delta_x,delta_y,visc)); % (eq 6.60, pag 330)
        end
    end

    % Zerar os devidos valores de F e G
    for j = 2:NJ-1
        % Paredes esquerda e direita
        F_mat(j,2) = 0; % (eq 6.94, pag 346)
        F_mat(j,NI) = 0;
    end
    for i = 2:NI-1
        % Paredes superior e inferior
        G_mat(2,i) = 0; % (eq 6.95, pag 346)
        G_mat(NJ,i) = 0;
    end

    % Loop para o cálculo das pressões
    for loop1 = 1:iter1
        disp(['PSOR - Iteração ',num2str(loop1)])

        % Guardar a matriz de dados antes das operações
        p_mat_old = p_mat;

        % Aplicar as condições de fronteira
        for j = 2:NJ-1
            % Paredes esquerda e direita
            p_mat(j,1) = p_mat(j,2); % (eq 6.97, pag 347)
            p_mat(j,NI) = p_mat(j,NI-1);
        end
        for i = 2:NI-1
            % Paredes superior e inferior
            p_mat(1,i) = p_mat(2,i); % (eq 6.98, pag 347)
            p_mat(NJ,i) = p_mat(NJ-1,i);
        end

        % Calcular os valores das pressões no interior da cavidade
        for i = 2:NI-1
            for j = 2:NJ-1
                % Paredes verticais
                if i == 2 || i == NI-1
                    val1 = u_mat(j,i+1);
                else
                    val1 = F_mat(j,i+1);
                end
                % Paredes horizontais
                if j == 2 || j == NJ-1
                    val2 = v_mat(j+1,i);
                else
                    val2 = G_mat(j+1,i);
                end
                % O motivo para fazer F = u e G = v se deve às equações (6.58)
                % e (6.59)
                fij = rho/delta_t*((val1-F_mat(j,i))/delta_x + (val2-G_mat(j,i))/delta_y);
%                if i == 2 || i == NI-1 || j == 2 || j == NJ-1
%                    fij = rho/delta_t*((u_mat(j,i+1)-F_mat(j,i))/delta_x + (v_mat(j+1,i)-G_mat(j,i))/delta_y);
%                else
%                    fij = rho/delta_t*((F_mat(j,i+1)-F_mat(j,i))/delta_x + (G_mat(j+1,i)-G_mat(j,i))/delta_y);
%                end
%                p_mat(j,i) = (1-w)*p_mat_old(j,i) + w/(2*(1+beta^2))*(p_mat_old(j,i+1) + p_mat_old(j,i-1) ...
%                             + beta^2*p_mat_old(j+1,i) + beta^2*p_mat_old(j-1,i) - delta_x^2*fij);
                p_mat(j,i) = (1-w)*p_mat(j,i) + w/(2*(1+beta^2))*(p_mat(j,i+1) + p_mat(j,i-1) ...
                             + beta^2*p_mat(j+1,i) + beta^2*p_mat(j-1,i) - delta_x^2*fij);
            end
        end

        % Normalizar as pressões
        p0 = p_mat(2,2);
        for i = 2:NI-1
            for j = 2:NJ-1
                p_mat(j,i) = p_mat(j,i) - p0;
            end
        end
        disp(p_mat),disp('')

        % Checar convergência
        R_mat = zeros(NJ,NI);
        for i = 2:NI-1
            for j = 2:NJ-1
                R_mat(j,i) = rho/delta_t*((F_mat(j,i+1) - F_mat(j,i))/delta_x + (G_mat(j+1,i) - G_mat(j,i))/delta_y) - ...
                             ((p_mat(j,i+1) - 2*p_mat(j,i) + p_mat(j,i-1))/delta_x^2 + (p_mat(j+1,i) - 2*p_mat(j,i) + p_mat(j-1,i))/delta_y^2);
            end
        end
        disp(norm(R_mat))
        if norm(R_mat) <= 1e-8
            disp('PSOR - |Convergência|')
            break
        end

    end

    % Guardar as velocidades para fins de comparação
    u_mat_old = u_mat;
    v_mat_old = v_mat;

    % Obter componentes u da velocidade
    for i = 2:NI-2
        for j = 2:NJ-1
            u_mat(j,i+1) = F_mat(j,i+1) - delta_t/rho*(p_mat(j,i+1) - p_mat(j,i))/delta_x; % (eq 6.103, pag 348)
        end
    end

    % Obter componentes v da velocidade
    for i = 2:NI-1
        for j = 2:NJ-2
            v_mat(j+1,i) = G_mat(j+1,i) - delta_t/rho*(p_mat(j+1,i) - p_mat(j,i))/delta_y; % (eq 6.104, pag 348)
        end
    end

    % Fazer o gráfico
    figure(1),clf
    quiver(X,Y,(u_mat(2:NJ-1,2:NI-1)),(v_mat(2:NJ-1,2:NI-1)),'r'),grid on
    title(['Iteração ',num2str(loop2)])

    % Checar a convergência da simulação completa
    S = 0;
    for i = 2:NI-1
        for j = 2:NJ-1
            S = S + abs(u_mat(j,i) - u_mat_old(j,i)) + abs(v_mat(j,i) - v_mat_old(j,i));
        end
    end
    if S/delta_t <= eps
        disp('Convergência total')
        break
    end

end

%% Obter resultantes
%res = zeros(NJ,NI);
%for i = 2:NI-1
%    for j = 2:NJ-1
%        res(j,i) = sqrt(u_mat(j,i)^2 + v_mat(j,i)^2);
%    end
%end
%
%% Adicionar ao gráfico
%hold on
%contourf(X,Y,res(2:NJ-1,2:NI-1));

% Fazer gráfico da pressão%figure(2),clf
%contourf(X,Y,p_mat(2:NJ-1,2:NI-1))


% Calcular função corrente
%psi = zeros(NJ,NI);
%for i = 2:NI-1
%    for j = 2:NJ-1
%        psi(j,i) = psi (j,i-1) + delta_y*u_mat(j,i);
%%        psi(j,i) = (u_mat(j,i) - u_mat(j-1,i))/delta_y - (v_mat() - v_mat())/delta_x;
%    end
%end
%figure(3),clf
%contour(X,Y,psi(2:NJ-1,2:NI-1)),grid on
