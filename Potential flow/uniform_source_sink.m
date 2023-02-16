clc,clear

% Obter componentes de velocidade
%syms u v X Y m pi
%psi = (u.*Y - v.*X) + (m/(2*pi)*atand(Y./X)); %(escoamento uniforme) + (fonte/sumidouro)
%disp(matlabFunction(diff(psi,Y))) % u
%disp(matlabFunction(diff(-psi,X))) % v

% Escoamento uniforme
alpha = 45;
V = 1;

% Fonte/sumidouro
m = 1;

num = 5;
vec_x = linspace(-10,10,num);
vec_y = linspace(-10,10,num);

[X,Y] = meshgrid(vec_x,vec_y);

u1 = V*cosd(alpha);
v1 = V*sind(alpha);
u = zeros(num);
v = zeros(num);
for i = 1:num
    for j = 1:num
        u(i,j) = u1 + 90 * m ./ (pi ^ 2 * X(i,j) .* (1 + Y(i,j) .^ 2 ./ X(i,j) .^ 2));
        v(i,j) = v1 + 90 * Y(i,j) .* m ./ (pi ^ 2 * X(i,j) .^ 2 .* (1 + Y(i,j) .^ 2 ./ X(i,j) .^ 2));
    end
end

figure(1),clf
quiver(X,Y,u,v,'r'),grid on,hold on
psi = (u.*Y - v.*X) + (m/(2*pi)*atand(Y./X));
contour(X,Y,psi,20)
