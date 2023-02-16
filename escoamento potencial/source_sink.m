clc,clear

% Fonte/sumidouro
m = 1;

num = 20;
vec_x = linspace(-10,10,num);
vec_y = linspace(-10,10,num);

[X,Y] = meshgrid(vec_x,vec_y);

u = zeros(num);
v = zeros(num);
for i = 1:num
    for j = 1:num
        u(i,j) = m*X(i,j)/(2*pi*(X(i,j)^2 + Y(i,j)^2));
        v(i,j) = m*Y(i,j)/(2*pi*(X(i,j)^2 + Y(i,j)^2));
    end
end

figure(1),clf
quiver(X,Y,u,v,'r'),grid on,hold on
%streamline(X,Y,u,v,X,Y)
psi =  m/(2*pi)*atand(Y./X);
contour(X,Y,psi,'b');

