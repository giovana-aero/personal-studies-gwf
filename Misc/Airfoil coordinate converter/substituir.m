function substituir(arquivo)

correction = fopen(arquivo,'r');
texto = fscanf(correction,'%c');
fclose(correction);

for i = 1:size(texto,2)
    if texto(i) == '.'
        texto(i) = ',';
    end
end

delete(arquivo);
correction = fopen(arquivo,'w');
fprintf(correction,'%c', texto);
fclose(correction);

end