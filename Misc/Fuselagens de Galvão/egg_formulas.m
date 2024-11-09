function egg_y = egg_formulas(x,B,L,n,op)

% x -> abscissas
% B -> maximum breath
% L -> egg length
% n -> número que altera o formato 
% op -> opção da fórmula a ser usada

w = (L-B)/2*n; % eqn 2


if op == 1
    egg_y = B/2*sqrt((L^2 - 4*x.^2)./(L^2 + 8*w*x + 4*w^2)); % eqn 1
elseif op == 2
    egg_y = B/2*sqrt(((L^2 - 4*x.^2)*L)./(2*(L - 2*w)*x.^2 + (L^2 + 8*L*w - 4*w^2)*x + 2*L*w^2 + L^2*w + L^3));


end