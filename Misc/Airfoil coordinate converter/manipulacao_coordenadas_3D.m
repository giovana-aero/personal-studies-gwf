clc,clear,fclose all;

% Manipulador de coordenadas de aerofólio
% Giovana Weffort
% Dezembro de 2022

% Dados ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Nota: o arquivo de coordenadas original deve estar formatado usando pontos
% como separadores decimais e não deve conter nenhum outro dado além das
% coordenadas em si (como o nome do perfil e o número de pontos)

% As rotações são feitas em termos da origem (isto é, o bordo de ataque).
% Devido a isso, tome o bordo de ataque como referência para as translações.

% Nome do arquivo (especificar também a extensão, como '.txt', ao final do nome)
nome = 'FX 61-163 AIRFOIL.txt';
% Reescala (alterar o comprimento da corda pra qual valor, sabendo que o
% original é unitário?)
scale = 1;
% Translações
dx = 0;
dy = 0;
dz = 0;
% Rotações (em graus)
thx = 0;
thy = 0;
thz = 0;

% Opções de impressão
formato = 0; % 1 para imprimir em sldcrv (outro valor para txt)
virgulaz = 0; % 1 para usar vírgulas como separador decimal
              % (outro valor para pontos)

% Leitura e conversão do formato das coordenadas ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Ler as coordenadas
coo = dlmread(nome);

% Modificar as coordenadas se necessário
% (o formato deve ser tal que o contorno comece no bordo de fuga, siga ao longo
% do extradorso, passe pelo bordo de ataque e volte ao bordo de fuga pelo
% intradorso)
delta = zeros(1,size(coo,1)-1);
for i = 1:(size(coo,1)-1)
    delta(i) = coo(i+1,1) - coo(i,1);
end

delta_b = delta >= 0;
sum_n = sum(delta_b==0);
if sum_n < size(coo,1)*0.1
    % Coordenadas em duas partes, ambas começando no bordo de ataque e
    % terminando no bordo de fuga
    for i = 1:size(delta,2)
        if delta_b(i) == 0
            zero_p = i;
            break
        end
    end

    coo = [flip(coo(2:zero_p,:));coo(zero_p+1:end,:)];
end

% Adicionar mais uma coluna pra completar as três dimensões
coo = [coo,zeros(size(coo,1),1)];

% Fazer uma cópia das coordenadas
coo_mod = coo;

% Manipulações das coordenadas ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Rotações
Rx = [1,0,0;
      0,cosd(thx),-sind(thx);
      0,sind(thx),cosd(thx)];
Ry = [cosd(thy),0,sind(thy);
      0,1,0;
      -sind(thy),0,cosd(thy)];
Rz = [cosd(thz),-sind(thz),0;
      sind(thz),cosd(thz),0;
      0,0,1];
R = Rz*Ry*Rx;

if ~isequal(eye(3),R)
    for i = 1:size(coo,1)
        coo_mod(i,:) = (R*coo(i,:)')';
    end
end

% Reescala
coo_mod = coo_mod*scale;

% Translações
coo_mod(:,1) = coo_mod(:,1) + dx;
coo_mod(:,2) = coo_mod(:,2) + dy;
coo_mod(:,3) = coo_mod(:,3) + dz;

% Visualizar a geometria e imprimir as coordenadas ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Fazer um gráfico
figure(1),clf
plot3(coo(:,1),coo(:,2),coo(:,3),'k'),grid on,hold on,axis equal
plot3(coo_mod(:,1),coo_mod(:,2),coo_mod(:,3),'r')
legend('Originais','Modificadas')
xlabel('x'),ylabel('y'),zlabel('z')

% Imprimir coordenadas
nome2 = nome(1:end-4);
if formato == 1
    nome2 = [nome2,'_mod.sldcrv'];
else
    nome2 = [nome2,'_mod.txt'];
end

fid = fopen(nome2,'w');
fprintf(fid,'%f %f %f\n',coo_mod');
fclose(fid);

if virgulaz == 1
    % Ler as coordenadas como uma string e substituir pontos por vírgulas
    fid = fopen(nome2,'r');
    text = fscanf(fid,'%c');
    text(text=='.') = ',';
    fclose(fid);

    % Sobrescrever o arquivo anterior
    fid = fopen(nome2,'w');
    fprintf(fid,'%s\n',text);
    fclose(fid);
end
