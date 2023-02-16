clc,clear

% Escoamento uniforme
alpha = 45;
V = 1;

num = 10;
vec_x = linspace(-10,10,num);
vec_y = linspace(-10,10,num);

[X,Y] = meshgrid(vec_x,vec_y);

u = zeros(num);
v = zeros(num);
for i = 1:num
    for j = 1:num
        u(i,j) = V*cosd(alpha);
        v(i,j) = V*sind(alpha);
    end
end
figure(1),clf
quiver(X,Y,u,v,'r'),grid on,hold on
%streamline(X,Y,u,v,X,Y);

psi = u.*Y - v.*X;
contour(X,Y,psi,'b')
