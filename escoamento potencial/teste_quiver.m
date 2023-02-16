clc,clear

spacing = 0.2;
[X,Y] = meshgrid(-2:spacing:2);

Z = X.*exp(-X.^2-Y.^2);
[DX,DY] = gradient(Z,spacing);

%figure(1),clf
%quiver(X,Y,DX,DY),hold on,grid on
%contour(X,Y,Z)



%psi = b*y - a*x % Corrente
%phi = b*x + a*y % Potencial

a = 1;
b = 1;

spacing = 0.1;
[X,Y] = meshgrid(-1:spacing:1);
% (exemplo 4.7 do white - mecânica dos fluidos)
u = 1*(X.^2 - Y.^2);
v = -2*1*X.*Y;
figure(1),clf
h1 = quiver(X,Y,u,v,'r');,grid on, hold on
%set(h1,'AutoScale','on','AutoScaleFactor',5)
[X2,Y2] = meshgrid(-1:0.2:1);
streamline(X,Y,u,v,X2,Y2);

