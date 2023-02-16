clc,clear

% Escoamento sobre degrau usando SOLA

% Dados do problema
h = 1; % [m]
H = 2*h; % [m]
g = 1.8; % [m]
G = 5; % [m]
% Água a 20 °C
visc = 0.000001; % [m^2/s] (viscosidade cinemática) https://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html
rho = 1e3; % [kg/m^3]
Re = 10;
u0 = Re*visc/h;

% Dados da malha
NI1 = 9;  % Horizontal: entrada até o final do degrau
NI2 = 16; % Horizontal: final do degrau até a saída
NJ1 = 5;  % Vertical: fronteira inferior até abaixo do degrau
NJ2 = 5;  % Vertical: acima do degrau até a fronteira superior
NI = NI1 + NI2;
NJ = NJ1 + NJ2;
%delta_x1 = g/NI1;
%delta_x2 = (G-g)/NI2;
%delta_y1 = h/NJ1;
%delta_y2 = (H-h)/NJ2;
delta_x = G/NI;
delta_y = H/(NJ-2);
% (por hora, pra facilitar, a malha será mantida uniforme e alterações nas
% dimensões da geometria não serão permitidas)

% Iterações
iter1 = 50; % Iterações das correções de pressão
iter2 = 10; % Iterações da solução completa
omega = 1.5; % Fator de relaxação

% Inicializar malhas
u_mat = zeros(NJ,NI);
v_mat = zeros(NJ,NI);
p_mat = zeros(NJ,NI);
deltau_mat = zeros(NJ,NI);
deltav_mat = zeros(NJ,NI);
deltap_mat = zeros(NJ,NI);
F_mat = zeros(NJ,NI);
G_mat = zeros(NJ,NI);
D_mat = zeros(NJ,NI);

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

% Loop principal
for loop2 = 1:iter2
    disp(['Solucionador - Iteração ',num2str(loop2)])

    % Atribuir as condições de fronteira
    % Velocidades na entrada
    u_mat(NJ1+1:NJ-1,1) = repmat(u0,size(u_mat(NJ1+1:NJ-1,1),1),1);
    % Condição de não escorregamento na parede de cima
    for i = 2:NI-1
        u_mat(NJ,i) = -u_mat(NJ-1,i);
    end
    % Condição de não escorregamento na parede de baixo (antes do degrau)
    for i = 2:NI1
        u_mat(NJ1,i) = u_mat(NJ1+1,i);
    end
    % Condição de não escorregamento na parede de baixo (depois do degrau)
    for i = NI+1:NI
        u_mat(1,i) = u_mat(2,i);
    end
    % Condição de não escorregamento na parede lateral do degrau
    for j = 2:NJ1
        v_mat(j,NI1) = v_mat(j,NI1+1);
    end
    % Velocidades normais nulas na parede de cima
    for i = 1:NI-1
        v_mat(NJ,i) = 0;
    end
    % Velocidades normais nulas na parede de baixo (antes do degrau)
    for i = 1:NI1
        v_mat(NJ1+1,i) = 0;
    end
    % Velocidades normais nulas na parede de baixo (depois do degrau)
    for i = NI1+1:NI-1
        v_mat(2,i) = 0;
    end
    % Velocidades normais nulas na parede lateral do degrau
    for j = 2:NJ1
        u_mat(j,NI+1) = 0;
    end

    % Obter delta_t
    delta_t = calc_delta_t(u_mat,v_mat,delta_x,delta_y,visc,tau);

    % Calcular F (Parte de cima)
    for i = 1:NI-1
        for j = NJ1+1:NJ-1
            F_mat(j,i+1) = u_mat(j,i+1) + delta_t*(-CONV_F(i,j,u_mat,v_mat,delta_x,delta_y,visc) + VISC_F(i,j,u_mat,delta_x,delta_y,visc)); % (eq 6.56, pag 329)
        end
    end

    % Calcular F (Parte de baixo - depois do degrau)
    for i = NI1+1:NI-2
        for j = 2:NJ1
            F_mat(j,i+1) = u_mat(j,i+1) + delta_t*(-CONV_F(i,j,u_mat,v_mat,delta_x,delta_y,visc) + VISC_F(i,j,u_mat,delta_x,delta_y,visc)); % (eq 6.56, pag 329)
        end
    end

    % Calcular G (Parte de cima)
    for i = 1:NI-1
        for j = NJ+1:NJ-2
            G_mat(j+1,i) = v_mat(j+1,i) + delta_t*(-CONV_G(i,j,u_mat,v_mat,delta_x,delta_y,visc) + VISC_G(i,j,v_mat,delta_x,delta_y,visc)); % (eq 6.60, pag 330)
        end
    end

    % Calcular G (Parte de baixo)
    for i = NI1+1:NI-1
        for j = 2:NJ1
            G_mat(j+1,i) = v_mat(j+1,i) + delta_t*(-CONV_G(i,j,u_mat,v_mat,delta_x,delta_y,visc) + VISC_G(i,j,v_mat,delta_x,delta_y,visc)); % (eq 6.60, pag 330)
        end
    end

    % Encontrar as velocidades no interior do canal
    % Componentes horizontais (parte de cima)
    for i = 1:NI-1
        for j = NJ1+1:NJ-1
            u_mat(j,i+1) = F_mat(j,i+1) - delta_t/rho*(p_mat(j,i+1) - p_mat(j,i))/delta_x;
        end
    end
    % Componentes horizontais (parte de baixo)
    for i = NI1+1:NI-2
        for j = 2:NJ1
            u_mat(j,i+1) = F_mat(j,i+1) - delta_t/rho*(p_mat(j,i+1) - p_mat(j,i))/delta_x;
        end
    end
    % Componentes verticais (parte de cima)
    for i = 1:NI-1
        for j = NJ+1:NJ-2
            v_mat(j+1,i) = G_mat(j+1,i) - delta_t/rho*(p_mat(j+1,i) - p_mat(j,i))/delta_y;
        end
    end
    % Componentes verticais (parte de baixo)
    for i = NI1+1:NI-1
        for j = 2:NJ1
            v_mat(j+1,i) = G_mat(j+1,i) - delta_t/rho*(p_mat(j+1,i) - p_mat(j,i))/delta_y;
        end
    end

    % Solução da pressão e mais correções das velocidades
    for loop1 = 1:iter1
        disp(['Correções - Iteração ',num2str(loop1)])

        % Calcular D (parte de cima)
        for i = 1:NI-1
            for j = NJ1+1:NJ-1
                D_mat(j,i) = ((u_mat(j,i+1) + deltau_mat(j,i+1)) - (u_mat(j,i) + deltau_mat(j,i)))/delta_x ...
                         + ((v_mat(j+1,i) + deltav_mat(j+1,i)) - (v_mat(j,i) + deltav_mat(j,i)))/delta_y;
            end
        end
        % Calcular D (parte de baixo )
        for i = NI1+1:NI-1
            for j = 2:NJ-1
                D_mat(j,i) = ((u_mat(j,i+1) + deltau_mat(j,i+1)) - (u_mat(j,i) + deltau_mat(j,i)))/delta_x ...
                         + ((v_mat(j+1,i) + deltav_mat(j+1,i)) - (v_mat(j,i) + deltav_mat(j,i)))/delta_y;
            end
        end

        % Calcular as correções de pressão
        % Parte de cima
        for i = 1:NI-1
            for j = NJ1+1:NJ-1
                deltap_mat(j,i) = -omega*D_mat(j,i)/(2*delta_t/rho*(delta_x^-2 + delta_y^-2));
            end
        end
        % Parte de baixo
        for i = NI1+1:NI-1
            for j = 2:NJ-1
                deltap_mat(j,i) = -omega*D_mat(j,i)/(2*delta_t/rho*(delta_x^-2 + delta_y^-2));
            end
        end

        % Corrigir velocidades
        % Parte de cima
        for i = 1:NI-1
            for j = NJ1+1:NJ-1
                u_mat(j,i+1) = u_mat(j,i+1) + delta_t/rho*deltap_mat(j,i)/delta_x;
                u_mat(j,i) = u_mat(j,i) - delta_t/rho*deltap_mat(j,i)/delta_x;
                v_mat(j+1,i) = v_mat(j+1,i) + delta_t/rho*deltap_mat(j,i)/delta_y;
                v_mat(j,i) = v_mat(j,i) - delta_t/rho*deltap_mat(j,i)/delta_y;
            end
        end
        % Parte de baixo
        for i = NI1+1:NI-1
            for j = 2:NJ-1
                u_mat(j,i+1) = u_mat(j,i+1) + delta_t/rho*deltap_mat(j,i)/delta_x;
                u_mat(j,i) = u_mat(j,i) - delta_t/rho*deltap_mat(j,i)/delta_x;
                v_mat(j+1,i) = v_mat(j+1,i) + delta_t/rho*deltap_mat(j,i)/delta_y;
                v_mat(j,i) = v_mat(j,i) - delta_t/rho*deltap_mat(j,i)/delta_y;
            end
        end

        disp(max(D_mat))



    end

end


