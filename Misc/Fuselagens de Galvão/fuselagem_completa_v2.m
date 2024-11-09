%%

% Giovana Weffort
% Dragonfly
% Fevereiro de 2022

% Esta versão do código trabalha gerando a fuselagem toda em um só script, ao
% invés de usar uma função pra cada 'metade' da fuselagem


%%
clear,clc,fclose('all');

%% Dados da parte da frente (fuselagem de Galvão) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Nomes dos arquivos de coordenadas
%nome_L = '64-206.txt';
%nome_CB = '64-208.txt';
nome_L = 'n64008a.dat';
nome_CB = 'n64108.dat';

% Imprimir coordenadas?
print = 0; % Pra imprimir: print = 1

% Fatores de reescala (comprimento)
size_f_L = 6400/1.4;%/1.8;
%size_f_CB = 6400;%*0.9;
size_f_CB = 6400*0.9;

% Fatores de reescala (espessura)
th_f_L = 1;
th_f_CB_ex = 1;
th_f_CB_in = 1;

% Inverter eixo x?
inv_x = 1; % Pra inverter: inv_x = 1

% Número de seções
n_sec = 15;

% Ponto ao longo do eixo x onde as seções acabam
p_x = 3.5556e3/2;
%p_x = 4e3;

% Pontos na seção transversal
n_elipse = 20;



%% Dados da parte de trás da fuselagem ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% (separada em duas partes: uma parte de transição das elipses à cauda, e a
% cauda em si, com seção transversal uniforme (circular ou elíptica))
p_x2 = p_x + 2000; % Final da transição pra cauda
p_x3 = p_x2 + 3000; % Final da cauda

dz1 = 'E'; % Posição vertical do início do centro do trecho de transição (opção E encaixa exatamente na
           % mesma posição que a última seção da parte da frente)

% Número de seções transversais do trecho de transição
n_sec2 = 15;

% Número de seções transversais da cauda
n_sec3 = 15;

% Raios do início da cauda/final do trecho de transição
ry_c2 = 100;
rz_c2 = 100;
dz2 = -00; % Posição vertical do centro do início da cauda

% Raios do final da cauda
ry_c3 = 50;
rz_c3 = 50;
dz3 = -00; % Posição vertical do centro do final da cauda

% Configuração do cosspace2
CP_factor1 = 0; % lados
CP_factor2 = 0; % cima e baixo



%% Cálculos referentes à parte da frente ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Lados
% Ler coordenadas do perfil dos lados
[coo_L,zero_p_L] = converter_function_xfoil(nome_L);

% Pegar coordenadas do extradorso
coo_ex_L = coo_L(1:zero_p_L,:);

% Multiplicar por 100 pra dar certo com o método de Galvão
coo_ex_L = coo_ex_L*100;

% Fazer spline
xx = linspace(0,100,100);
spl_L_half = spline(coo_ex_L(1:zero_p_L,1),coo_ex_L(1:zero_p_L,2),xx);

% Transformar em fuselagem de Galvão
gv_L_half = zeros(1,size(spl_L_half,2));
for i = 1:size(spl_L_half,2)
    if spl_L_half(i) < 0
        gv_L_half(i) = -(-spl_L_half(i))^(3/2);
    else
        gv_L_half(i) = spl_L_half(i)^(3/2);
    end
end

% Aplicar reescala de espessura
gv_L_half = gv_L_half*th_f_L;

% Gerar matriz de coordenadas da fuselagem (lados)
coo_gv_L = [flip(xx') flip(gv_L_half');
    xx(:,2:end)' -gv_L_half(:,2:end)'];

% Deixar com comprimento unitário
coo_gv_L = coo_gv_L/100;

%% Cima/baixo
% Ler coordenadas do perfil de cima/baixo
[coo_CB,zero_p_CB] = converter_function_xfoil(nome_CB);

% Pegar coordenadas do extradorso e intradorso
coo_ex_CB = coo_CB(1:zero_p_CB,:);
coo_in_CB = coo_CB(zero_p_CB:end,:);

% Multiplicar por 100 pra dar certo com o método de Galvão
coo_ex_CB = coo_ex_CB*100;
coo_in_CB = coo_in_CB*100;

% Fazer spline
xx = linspace(0,100,100);
spl_CB_ex = spline(coo_ex_CB(1:zero_p_CB,1),coo_ex_CB(1:zero_p_CB,2),xx);
spl_CB_in = spline(coo_in_CB(1:zero_p_CB,1),coo_in_CB(1:zero_p_CB,2),xx);

% Transformar em fuselagem de galvão
gv_CB_ex = zeros(1,size(spl_CB_ex,2));
gv_CB_in = zeros(1,size(spl_CB_in,2));
for i = 1:size(spl_CB_ex,2)
    if spl_CB_ex(i) < 0
        gv_CB_ex(i) = -(-spl_CB_ex(i))^(3/2);
    else
        gv_CB_ex(i) = spl_CB_ex(i)^(3/2);
    end
end
for i = 1:size(spl_CB_in,2)
    if spl_CB_in(i) < 0
        gv_CB_in(i) = -(-spl_CB_in(i))^(3/2);
    else
        gv_CB_in(i) = spl_CB_in(i)^(3/2);
    end
end

% Aplicar reescala de espessura
gv_CB_ex = gv_CB_ex*th_f_CB_ex;
gv_CB_in = gv_CB_in*th_f_CB_in;

coo_gv_CB = [flip(xx') flip(gv_CB_ex');
            xx(:,2:end)' gv_CB_in(:,2:end)'];

% Achar a linha de curvatura de cima/baixo
camber = (gv_CB_ex + gv_CB_in)/2;

% Deixar com comprimento unitário
coo_gv_CB = coo_gv_CB/100;

%% Configurar os planos de cada contorno
camber_pl = [camber flip(camber(1:end-1))]/100;

% Perfil da lateral fica no plano XY
coo_gv_L_pl = [coo_gv_L camber_pl'];

% Perfil de cima/baixo fica no plano XZ
coo_gv_CB_pl = [coo_gv_CB(:,1) zeros(size(coo_gv_CB,1),1) coo_gv_CB(:,2)];
% coo_gv_CB_pl = [coo_gv_CB(:,1) camber_pl' coo_gv_CB(:,2)];

%% Gráfico 3D dos contornos (sem modificações)
%figure(1),clf,hold on,axis equal,grid on
%plot3(coo_gv_L_pl(:,1),coo_gv_L_pl(:,2),coo_gv_L_pl(:,3))
%plot3(coo_gv_CB_pl(:,1),coo_gv_CB_pl(:,2),coo_gv_CB_pl(:,3))
%title('fuselagem sem modificaÃ§Ãµes')

%% Gráficos pra checar se está tudo em ordem
%figure(10),clf
%plot(coo_gv_L(:,1),coo_gv_L(:,2))
%grid on,axis equal
%
%figure(11),clf
%plot(coo_gv_CB(:,1),coo_gv_CB(:,2))
%grid on,axis equal


%% Fazer seçõees transversais no início da fuselagem
coo_gv_L_pl_mod = [coo_gv_L_pl(:,1)*-1 coo_gv_L_pl(:,2:3)]*size_f_L;
coo_gv_CB_pl_mod = [coo_gv_CB_pl(:,1)*-1 coo_gv_CB_pl(:,2:3)]*size_f_CB;

% Gráfico 3D dos contornos (modificados)
figure(20),clf,grid on,hold on,axis equal
%plot3(coo_gv_L_pl_mod(:,1),coo_gv_L_pl_mod(:,2),coo_gv_L_pl_mod(:,3))
%plot3(coo_gv_CB_pl_mod(:,1),coo_gv_CB_pl_mod(:,2),coo_gv_CB_pl_mod(:,3))
title('fuselagem modificada')
xlabel('x'),ylabel('y'),zlabel('z')

% Criar splines
if inv_x == 1
    p_x = -p_x;
    xx = -xx;
    xx2 = linspace(0,p_x,n_sec); % distribuição de seções

    % Fazer splines que vão até o ponto de interesse (p_x)
    spl_gv_L_half = spline(xx/100*size_f_L,gv_L_half/100,xx2);
    spl_gv_CB_ex = spline(xx/100*size_f_CB,gv_CB_ex/100,xx2);
    spl_gv_CB_in = spline(xx/100*size_f_CB,gv_CB_in/100,xx2);
    spl_camber = spline(xx/100*size_f_CB,camber/100,xx2);

end

% Fazer um gráfico da fuselagem cortada pra checar
% scatter3(xx2,spl_gv_L_half,zeros(1,length(xx2)))
% scatter3(xx2,-spl_gv_L_half,zeros(1,length(xx2)))
% scatter3(xx2,zeros(1,length(xx2)),spl_gv_CB_ex)
% scatter3(xx2,zeros(1,length(xx2)),spl_gv_CB_in)
%plot3(xx/100*size_f_CB,zeros(1,length(camber)),camber/100*size_f_CB,'y--')

% Vetor com informações das elipses
delta_z = spl_gv_CB_ex-spl_camber;

% Struct pra guardar os pontos das elipses
empty.coo = []; %zeros(n_elipse,3);
elipsez = repmat(empty,n_sec,1);

% Criar elipses (plano YZ)
for i = 2:n_sec

     %posição do centro da elipse
     yCenter = 0;
     zCenter = spl_camber(i);

     %raios da elipse
     yRadius = spl_gv_L_half(i)*size_f_L;
     zRadius = delta_z(i)*size_f_CB;

     %geração da elipse baseando-se em um número n de pontos
     theta = linspace(0,360,n_elipse+1)+90;
     theta = flip(theta);

     y = yRadius * cosd(theta) + yCenter;
 %     y(:,n_elipse+1)=[]; %removendo a última linha dessa matriz
     z = zRadius * sind(theta) + zCenter;
 %     z(:,n_elipse+1)=[]; %removendo a última linha dessa matriz

     elipsez(i).coo = [repelem(xx2(i),n_elipse+1)' y' z'+spl_camber(i)*size_f_CB];

end

%% Visualizar as elipses
%figure(1),clf
for i = 2:n_sec
    scatter3(elipsez(i).coo(:,1),elipsez(i).coo(:,2),elipsez(i).coo(:,3))
end
hold on,grid on,axis equal,xlabel('x'),ylabel('y'),zlabel('z')


%% Cálculos referentes à parte de trás ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% << Trecho de transição >>
% Recalcular a posição vertical do início do trecho
if dz1 == 'E'
    dz1 = spl_camber(n_sec);
end

% Distribuições relevantes
if inv_x == 1
    x2 = linspace(p_x,-p_x2,n_sec2);
end
v_pos2 = linspace(dz1+spl_camber(n_sec)*size_f_CB,dz2,n_sec2); % Posições verticais do centro
%ry_range2 = linspace(spl_gv_L_half(n_sec)*size_f_L,ry_c2,n_sec2); % Raios horizontais
ry_range2 = cosspace2(spl_gv_L_half(n_sec)*size_f_L,ry_c2,n_sec2,CP_factor1);
%rz_range2 = linspace(delta_z(n_sec)*size_f_CB,rz_c2,n_sec2);
rz_range2 = cosspace2(delta_z(n_sec)*size_f_CB,rz_c2,n_sec2,CP_factor2);

% Struct pra guardar os pontos das elipses
%empty.coo = []; %zeros(n_elipse,3);
elipsez2 = repmat(empty,n_sec2,1);

% Criar elipses (plano YZ)
for i = 2:n_sec2

     %posição do centro da elipse
     yCenter = 0;
     zCenter = v_pos2(i);

     %raios da elipse
     yRadius = ry_range2(i);
     zRadius = rz_range2(i);

     %geração da elipse baseando-se em um número n de pontos
     theta = linspace(0,360,n_elipse+1)+90;
     theta = flip(theta);

     y = yRadius * cosd(theta) + yCenter;
 %     y(:,n_elipse+1)=[]; %removendo a última linha dessa matriz
     z = zRadius * sind(theta) + zCenter;
 %     z(:,n_elipse+1)=[]; %removendo a última linha dessa matriz

     elipsez2(i).coo = [repelem(x2(i),n_elipse+1)',y',z'];

end

% Visualizar a seção de transição
for i = 2:n_sec2
    scatter3(elipsez2(i).coo(:,1),elipsez2(i).coo(:,2),elipsez2(i).coo(:,3))
end


% << Cauda >>
% Recalcular a posição vertical do início do trecho
%if dz2 == 'E'
%    dz2 = v_pos2(end);
%end

% Distribuições relevantes
if inv_x == 1
    x3 = linspace(-p_x2,-p_x3,n_sec3);
end
v_pos3 = linspace(dz2,dz3,n_sec3); % Posições verticais do centro
ry_range3 = linspace(ry_c2,ry_c3,n_sec3); % Raios horizontais
rz_range3 = linspace(rz_c2,rz_c3,n_sec3);

% Struct pra guardar os pontos das elipses
%empty.coo = []; %zeros(n_elipse,3);
elipsez3 = repmat(empty,n_sec3,1);

% Criar elipses (plano YZ)
for i = 2:n_sec3

     %posição do centro da elipse
     yCenter = 0;
     zCenter = v_pos3(i);

     %raios da elipse
     yRadius = ry_range3(i);
     zRadius = rz_range3(i);

     %geração da elipse baseando-se em um número n de pontos
     theta = linspace(0,360,n_elipse+1)+90;
     theta = flip(theta);

     y = yRadius * cosd(theta) + yCenter;
 %     y(:,n_elipse+1)=[]; %removendo a última linha dessa matriz
     z = zRadius * sind(theta) + zCenter;
 %     z(:,n_elipse+1)=[]; %removendo a última linha dessa matriz

     elipsez3(i).coo = [repelem(x3(i),n_elipse+1)',y',z'];

end

% Visualizar a seção de transição
for i = 2:n_sec3
    scatter3(elipsez3(i).coo(:,1),elipsez3(i).coo(:,2),elipsez3(i).coo(:,3))
end




%% Imprimir arquivo de coordenadas ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if print == 1

%    % Cabeçalho: indicar nomes dos aerofólios e data/hora na primeira linha
%    % do arquivo. O primeiro nome é do aerofólio que faz os contornos
%    % laterais da fuselagem e o segundo nome é o que faz os contornos de
%    % cima e de baixo
%    header = [nome_L(1:end-4) '|' nome_CB(1:end-4) '|' datestr(now)];
%
%    % Imprimir tudo
%    coo_gv_L_pl_mod = [coo_gv_L_pl(:,1)*-1 coo_gv_L_pl(:,2:3)];
%    coo_gv_CB_pl_mod = [coo_gv_CB_pl(:,1)*-1 coo_gv_CB_pl(:,2:3)];
%
%    txt_nome = 'coordenadas_fuselagem_metre.txt';
%    arquivo = fopen(txt_nome,'w');
%    fprintf(arquivo,'%s\n',header);
%
%    % Imprimir perfis
%    fprintf(arquivo,'%f %f %f\n', size_f_L*coo_gv_L_pl_mod');
%    fprintf(arquivo,'~~~~~~~~~~~~~~~~~~~~\n');
%    fprintf(arquivo,'%f %f %f\n', size_f_CB*coo_gv_CB_pl_mod');

    % É necesário que não hajam caracteres especiais (acentos e afins)
    place = 'C:\Users\Guga Weffort\Desktop\fuselagens\';
     % Imprimir elipses (frente)
     for i = 1:n_sec
%         fprintf(arquivo,'EndCurve\nStartCurve\n');
         name = [place,'part1_section',num2str(i),'.sldcrv'];
         arquivo = fopen(name,'w');
         fprintf(arquivo,'%f %f %f\n',elipsez(i).coo');
         fclose(arquivo);
     end
     % Imprimir elipses (transição)
     for i = 1:n_sec2
%         fprintf(arquivo,'EndCurve\nStartCurve\n');
         name = [place,'part2_section',num2str(i),'.sldcrv'];
         arquivo = fopen(name,'w');
         fprintf(arquivo,'%f %f %f\n',elipsez2(i).coo');
         fclose(arquivo);
     end
     % Imprimir elipses (cauda)
     for i = 1:n_sec3
%         fprintf(arquivo,'EndCurve\nStartCurve\n');
         name = [place,'part3_section',num2str(i),'.sldcrv'];
         arquivo = fopen(name,'w');
         fprintf(arquivo,'%f %f %f\n',elipsez3(i).coo');
         fclose(arquivo);
     end


%    fclose(arquivo);

end
