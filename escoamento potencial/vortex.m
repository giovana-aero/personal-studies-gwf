clc,clear

%% Obter dados da velocidade
%syms gamma pi r r0 theta
%psi = -gamma/(2*pi)*log(r/r0);
%disp(matlabFunction(diff(1/r*psi,theta))) % v_r
%disp(matlabFunction(diff(-psi,r))) % v_theta

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

    end
end
