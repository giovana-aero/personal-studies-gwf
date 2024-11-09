function plot_airfoil(coordenadas,edit,zero_p,nome)

switch edit
    case {1,2}
        plot(coordenadas(1:zero_p,1),coordenadas(1:zero_p,2),'k'),grid on,hold on
        plot(coordenadas(zero_p+1:end,1),coordenadas(zero_p+1:end,2),'k'),grid on,hold on
        axis square, set(gca,'ylim',[-.5 .5])
        title(nome(1:end-4))
        
    case 3
        plot(coordenadas(:,1),coordenadas(:,2),'k'),grid on
        axis square, set(gca,'ylim',[-.5 .5])
        title(nome(1:end-4))
        
end


end