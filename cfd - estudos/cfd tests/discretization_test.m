clc,clear

% Funções necessárias pras interpolações
function ans = sinal(x)
    if x >= 0
        ans = 1;
    else
        ans = -1;
    end
end
function ans = calc_pe(ui,delta_x,visc)
    ans = ui*delta_x/visc;
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

% Distribuições quaisquer de velocidades, representando valores na malha
u_mat = [0.525929   0.973689   0.092271   0.682930   0.450989;
         0.775562   0.146131   0.792923   0.116863   0.280581;
         0.842212   0.561598   0.518301   0.673664   0.590758;
         0.658519   0.717314   0.704562   0.449992   0.634621;
         0.798848   0.267554   0.865009   0.371052   0.892896];
v_mat = [0.652789   0.061395   0.773436   0.464222   0.346961;
         0.988365   0.137554   0.272983   0.039997   0.581730;
         0.851676   0.067959   0.553466   0.880488   0.753248;
         0.380543   0.755928   0.015876   0.661146   0.275922;
         0.275992   0.910133   0.877527   0.476625   0.258105];

% Dados
delta_x = 0.1;
delta_y = 0.1;
visc = 0.000001;

% Ponto na malha
i = 2;
j = 2;

% Inicializar matrizes dos termos
%conv_F = zeros(5);
%visc_F = zeros(5);
%conv_G = zeros(5);
%visc_G = zeros(5);

% << Direção x >>

% Derivadas dos termos viscosos
%VISC = (u_mat(j,i) - 2*u_mat(j,i+1) + u_mat(j,i+2))/delta_x^2;
%VISC = visc*(u_mat(j,i) - 2*u_mat(j,i+1) + u_mat(j,i+2))/delta_x^2 + visc*(u_mat(j+1,i+1) - 2*u_mat(j,i+1) + u_mat(j-1,i+1))/delta_y^2;
%disp(VISC)



% Derivadas dos termos convectivos

% del(u^2)/(del x)
%u_bar -> Velocidade de convecção (média aritmética)
%u_prp -> Propriedade transportada (interpolação)

% Velocidades de convecção
u_bar_ip3o2 = (u_mat(j,i+1) + u_mat(j,i+2))/2;
u_bar_ip1o2 = (u_mat(j,i) + u_mat(j,i+1))/2;
% Interpolações
ue = (u_mat(j,i+1) + u_mat(j,i+2))/2; % Velocidade estimada
pe = calc_pe(ue,delta_x,visc); % Calcular número de péclet
u_prp_ip3o2 = hybrid_itp(pe,ue,u_mat(j,i+1),u_mat(j,i+2)); % Interpolar

ue = (u_mat(j,i) + u_mat(j,i+1))/2; % Velocidade estimada
pe = calc_pe(ue,delta_x,visc); % Calcular número de péclet
u_prp_ip1o2 = hybrid_itp(pe,ue,u_mat(j,i),u_mat(j,i+1)); % Interpolar

ans = (u_bar_ip3o2*u_prp_ip3o2 - u_bar_ip1o2*u_prp_ip1o2)/delta_x;
disp(ans)

