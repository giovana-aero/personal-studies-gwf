clf,clear,clc,fclose('all');

% implementar fechamento de bordo de fuga no catia também


% testar a mudança de eixo com as outras funções


%% 

% Giovana Weffort 
% Dragonfly
% Março de 2021


%% Cabeçalho e opções
% Nota: o arquivo original deve conter apenas as coordenadas para que não
% hajam erros.

% Padrões de marcadores de casas decimais no meu pc:
% XFOIL: pontos
% ANSYS: vírgulas
% CATIA: pontos

% Opções
% 1 - XFOIL, XFLR5 e OpenVSP
% 2 - ANSYS
% 3 - CATIA e SolidWorks
op = 1;

% [ANSYS] Gerar uma única linha ou gerar extradorso e intradorso
% separadamente? [0/1]
contorno = 1;

% [ANSYS] Para perfis abertos: fechar o bordo de fuga? [0/1]
close_TE = 1;

% [ANSYS e CATIA] Trocar pontos por vírgulas? [0/1]
virgulaz = 1;

% [ANSYS e CATIA] Fator de reescala:
size_f = 1;

% [ANSYS e CATIA] Trocar o plano onde o perfil é gerado?
% 1 - xy
% 2 - yz
% 3 - zx
plano = 1;
inv_axis = 0; % Inverter eixos [0/1]

% [ANSYS e CATIA] Distância do plano
d_plano = 0;

% [ANSYS e CATIA] Inverter (espelhar) perfis? [0/1]
inv_horizontal = 0;
inv_vertical = 0;

% [ANSYS e CATIA] Mover o perfil verticalmente e/ou horizontalmente? [0/1]
d_horizontal = 0;
d_vertical = 0;


%%

% Ler coordenadas do arquivo texto
% nome = 'NACA 63012A.txt';
% nome = '64-012.txt';
% nome = 'RAF-48.txt';
% nome = 'FX60-100.txt';
% nome = 'e818.txt';
% nome = '2410.txt';
% nome = 'rutan wing airfoil.txt';
% nome = '0012.txt';
% nome = 'Wortmann FX 61-184.dat';
% nome = 'NACA 5 DIGITS.txt';
%nome = '5311.txt';
nome = 'FX 61-163 AIRFOIL.txt';
nome_perfil = nome(1:end-4);
coordenadas = dlmread(nome);

delta = zeros(1,size(coordenadas,1)-1);
for i = 1:(size(coordenadas,1)-1)
    delta(i) = coordenadas(i+1,1) - coordenadas(i,1);
end

delta_b = delta >= 0;
sum_p = sum(delta_b==1);
sum_n = sum(delta_b==0);


if sum_p < size(coordenadas,1)*0.1
    % Coordenadas em duas partes, ambas começando no bordo de fuga e
    % terminando no bordo de ataque
    edit = 1;
    for i = 1:size(delta,2)
        if delta_b(i) == 0
            zero_p = i;
            break
        end
    end
elseif sum_n < size(coordenadas,1)*0.1
    % Coordenadas em duas partes, ambas começando no bordo de ataque e
    % terminando no bordo de fuga
    edit = 2;
    for i = 1:size(delta,2)
        if delta_b(i) == 0
            zero_p = i;
            break
        end
    end
else
    % Coordenadas em uma única parte (geralmente começando e terminando no
    % bordo de fuga
    edit = 3; 
    
    if delta_b(1) == 0
        for i = 1:size(delta_b,2)
            if delta_b(i) == 1
                zero_p = i;
            break
            end
        end
        
    else
        for i = 1:size(delta_b,2)
            if delta_b(i) == 0
                zero_p = i;
                break
            end
        end
    end
    
end

plot_airfoil(coordenadas,edit,zero_p,nome)

%% Conversão das coordenadas
switch op
    case 1 % XFOIL
        
        arquivo = strcat(nome_perfil,'_XFOIL.dat');
        out = fopen(arquivo,'w');
        
        switch edit
            case 1
                
            case 2
                % Inverter a ordem da primeira parte
                coordenadas = [flip(coordenadas(2:zero_p,:)); coordenadas(zero_p+1:end,:)];
                
            case 3
                if coordenadas(1,1) == 1 
                    % Não é necessário editar as coordenadas

                end
            
        end
        
        fprintf(out,'%s\n',nome_perfil);
        fprintf(out,'%f %f\n',coordenadas');
        fclose(out);
        
        
        
    case 2 % ANSYS
        arquivo = strcat(nome_perfil,'_ANSYS.txt');
        out = fopen(arquivo,'w');
        v = 'n';
        
        switch edit
            case 1
                
            case 2
                if contorno == 1
                    % Inverter a ordem da primeira parte
                    coordenadas = [flip(coordenadas(2:zero_p,:)); coordenadas(zero_p+1:end,:)];
                    
                    % Se as coordenadas não geram um contorno fechado
                    if coordenadas(1,2) ~= coordenadas(end,2)
                        n1 = zeros(size(coordenadas,1),1) + 1;
                        n2 = (1:(size(coordenadas,1)))';
%                         n3 = zeros(size(coordenadas,1),1);
                        
                        if close_TE == 1
                            v = [2 1 coordenadas(1,:) 0; 2 2 coordenadas(end,:) 0];
                        end
                        
                    % Se as coordenadas geram um contorno fechado
                    else
                        n1 = zeros(size(coordenadas,1),1) + 1;
                        n2 = 1:(size(coordenadas,1)-1); n2 = [n2';0];
%                         n3 = zeros(size(coordenadas,1),1);
                        
                    end
                    
                else
                    n1 = [(zeros(size(coordenadas,1)/2,1)+1); (zeros(size(coordenadas,1)/2,1)+2)];
                    n2 = [1:(size(coordenadas,1)/2) 1:(size(coordenadas,1)/2)]';
%                     n3 = zeros(size(coordenadas,1),1);
                    
                    % Caso o bordo de fuga esteja aberto
                    if coordenadas(zero_p,2) ~= coordenadas(end,2)
                       if close_TE == 1
                           v = [3 1 coordenadas(zero_p,:) 0; 3 2 coordenadas(end,:) 0];
                       end
                    end
                    
                end
                
            case 3
                if coordenadas(1,1) == 1 
                    % Não é necessário editar as coordenadas
                    if contorno == 1
                        % Se as coordenadas não geram um contorno fechado
                        if coordenadas(1,2) ~= coordenadas(end,2)
                            n1 = zeros(size(coordenadas,1),1) + 1;
                            n2 = (1:(size(coordenadas,1)))';
%                             n3 = zeros(size(coordenadas,1),1);
                            
                        % Se as coordenadas geram um contorno fechado
                        else
                            n1 = zeros(size(coordenadas,1),1) + 1;
                            n2 = 1:(size(coordenadas,1)-1); n2 = [n2';0];
%                             n3 = zeros(size(coordenadas,1),1);
                            
                        end
                    else
                        coordenadas = [coordenadas(1:zero_p,:); coordenadas(zero_p:end,:)];
                        n1 = [(zeros(size(coordenadas,1)/2,1)+1); (zeros(size(coordenadas,1)/2,1)+2)];
                        n2 = [1:(size(coordenadas,1)/2) 1:(size(coordenadas,1)/2)]';
%                         n3 = zeros(size(coordenadas,1),1);
                        
                    end

                end
            
        end
        
        % Aplicar reescala
        coordenadas = coordenadas*size_f;
        
        % Aplicar inversões
        if inv_horizontal == 1
            coordenadas(:,2) = coordenadas(:,2)*(-1);
        end
        if inv_vertical == 1
            coordenadas(:,1) = coordenadas(:,1)*(-1);
        end
        
        % Aplicar deslocamentos horizontais/verticais
        coordenadas(:,1) = coordenadas(:,1) + d_horizontal;
        coordenadas(:,2) = coordenadas(:,2) + d_vertical;
        
        % Mudar plano e aplicar distância a ele
        switch plano
            case 1
                coordenadas = [coordenadas zeros(size(coordenadas,1),1)+d_plano];
            case 2
                coordenadas = [zeros(size(coordenadas,1),1)+d_plano coordenadas];
            case 3
                coordenadas = [coordenadas(:,2) zeros(size(coordenadas,1),1)+d_plano coordenadas(:,1)];
            otherwise
                coordenadas = [coordenadas zeros(size(coordenadas,1),1)+d_plano];
        end
        
        coordenadas = [n1 n2 coordenadas];
%         coordenadas = coordenadas*size_f;
        if isa(v,'double'),coordenadas = [coordenadas; v];end
        fprintf(out,'%d %d %f %f %f\n',coordenadas');
        fclose(out);
        
        
        
    case 3 % CATIA
        arquivo = strcat(nome_perfil,'_CATIA.txt');
        out = fopen(arquivo,'w');
        
        switch edit
            case 1
                
            case 2
                % Inverter a ordem da primeira parte
                coordenadas = [flip(coordenadas(2:zero_p,:)); coordenadas(zero_p+1:end,:)];
                
            case 3
                if coordenadas(1,1) == 1 
                    % Não é necessário editar as coordenadas

                end
            
        end
        
        % Aplicar reescala
        coordenadas = coordenadas*size_f;
        
        % Aplicar inversões
        if inv_horizontal == 1
            coordenadas(:,2) = coordenadas(:,2)*(-1);
        end
        if inv_vertical == 1
            coordenadas(:,1) = coordenadas(:,1)*(-1);
        end
        
        % Aplicar deslocamentos horizontais/verticais
        coordenadas(:,1) = coordenadas(:,1) + d_horizontal;
        coordenadas(:,2) = coordenadas(:,2) + d_vertical;
        
        % Mudar plano e aplicar distância a ele
        switch plano
            case 1
                coordenadas = [coordenadas zeros(size(coordenadas,1),1)+d_plano];
                if inv_axis == 1
                    coordenadas = [coordenadas(:,2) coordenadas(:,1) coordenadas(:,3)];
                end
            case 2
                coordenadas = [zeros(size(coordenadas,1),1)+d_plano coordenadas];
                if inv_axis == 1
                    coordenadas = [coordenadas(:,1) coordenadas(:,3) coordenadas(:,2)];
                end
            case 3
                coordenadas = [coordenadas(:,2) zeros(size(coordenadas,1),1)+d_plano coordenadas(:,1)];
                if inv_axis == 1
                    coordenadas = [coordenadas(:,3) coordenadas(:,2) coordenadas(:,1)];
                end
            otherwise
                coordenadas = [coordenadas zeros(size(coordenadas,1),1)+d_plano];
        end
        
        % Imprimir pro arquivo de texto
        fprintf(out,'%f %f %f\n',coordenadas');
        fclose(out);
        
    otherwise
        disp('Opção inválida')
        
end

% Trocar pontos por vírgulas (operação válida apenas pro ANSYS e CATIA)
if (virgulaz == 1 && op ~= 1), substituir(arquivo), end
    


