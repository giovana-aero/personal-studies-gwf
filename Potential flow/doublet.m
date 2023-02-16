clc,clear

% Obter dados da velocidade
%syms m s X Y psi pi
%psi = -(2*m*s)*Y/(2*pi*(X.^2 + Y.^2));
%disp(matlabFunction(diff(psi,Y))) % u
%disp(matlabFunction(diff(-psi,X))) % v

% Dipolo
m = 1;
s = 1;

num = 5;
vec_x = linspace(-10,10,num);
vec_y = linspace(-10,10,num);

[X,Y] = meshgrid(vec_x,vec_y);

u = zeros(num);
v = zeros(num);
for i = 1:num
    for j = 1:num
        u(i,j) = 2 * Y(i,j) .^ 2 .* m .* s ./ (pi * (X(i,j) .^ 2 + Y(i,j) .^ 2) .^ 2) - m .* s ./ (pi * (X(i,j) .^ 2 + Y(i,j) .^ 2));
        v(i,j) = -2 * X(i,j) .* Y(i,j) .* m .* s ./ (pi * (X(i,j) .^ 2 + Y(i,j) .^ 2) .^ 2);
    end
end

psi = -(2*m*s)*Y./(2*pi*(X.^2 + Y.^2));

figure(1),clf
quiver(X,Y,u,v,'r'),grid on,hold on
%contour(X,Y,psi,'b',50)
