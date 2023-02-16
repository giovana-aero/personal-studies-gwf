clc,clear

%(page 437)

% Domain data
L = 0.5; % [ft]
D = 0.01; % [ft]
ue = 1; % [ft/s]

% Fluid properties
rho = 0.002377; % [slug/ft^3]
mu = 3.737e-7; % [slug/ft/s]   (http://www.aerodynamics4students.com/properties-of-the-atmosphere/sea-level-conditions.php)
Re = ue*rho*D/mu;

% Solution configuration
delta_t = 0.001;
iter1 = 200; % Main loop
iter2 = 200; % Pressure estimations
tol1 = eps; % Convergence tolerance - main loop
tol2 = 1e-7; % Convergence tolerance - Pressure estimation
alpha_p = 0.1; % Underrelaxation factor

% Mesh configuration
NX = 21;
NY = 11;
delta_x = L/(NX-1);
delta_y = D/(NY-1);
NXu = NX+1;
NYu = NY;
NXv = NX+2;
NYv = NY+1;

% Initialize matrices
p_mat = zeros(NY,NX);
u_mat = zeros(NYu,NXu);
v_mat = zeros(NYv,NXv);
rhou = zeros(NYu,NXu);
rhov = zeros(NYv,NXv);

% Apply boundary condition at the top wall
u_mat(1,:) = u_mat(1,:) + ue;

% Arbitrary velocity value to test the quality of the solution
v_mat(5,15) = 0.5;
%u_mat(5,6) = -1;
%u_mat(2,16) = 10;
%v_mat(9,8) = 5;

% History of v at the above point
v_his = zeros(1,iter1);

for loop1 = 1:iter1

    disp(['Iteration ',num2str(loop1)])

    % Obtain rho*u
    for i = 1:NX-1
        for j = 2:NY-1
            v1 = (v_mat(j+1,i+1)+v_mat(j+1,i+2))/2;
            v2 = (v_mat(j,i+1)+v_mat(j,i+2))/2;
            A = -rho*((u_mat(j,i+2)^2-u_mat(j,i)^2)/(2*delta_x) + (u_mat(j+1,i+1)*v1-u_mat(j-1,i+1)*v2)/(2*delta_y)) +...
                mu*((u_mat(j,i+2)-2*u_mat(j,i+1)+u_mat(j,i))/delta_x^2 + (u_mat(j+1,i+1)-2*u_mat(j,i+1)+u_mat(j-1,i+1))/delta_y^2);

            rhou(j,i+1) = rhou(j,i+1) + A*delta_t - delta_t/delta_x*(p_mat(j,i+1)-p_mat(j,i));

            % Obtain velocities
            u_mat(j,i+1) = rhou(j,i+1)/rho;

        end
    end

    % Keep the value of the velocity at that one point
    v_his(loop1) = v_mat(5,15);

    % Changes to the discretization so the algorithm stays within bounds:
    % del^2(u)/del(x)^2 -> change to central differences (second order)
    % del(rho*u)/del(x) -> change to central differences (second order)
    % del(p)/del(x) -> change to rearward differences (first order)
%    i = NX;
%    for j = 2:NY-1
%        v1 = (v_mat(j+1,i+1)+v_mat(j+1,i+2))/2;
%        v2 = (v_mat(j,i+1)+v_mat(j,i+2))/2;
%        A = -rho*((u_mat(j,i+1)^2-u_mat(j,i-1)^2)/(2*delta_x) + (u_mat(j+1,i+1)*v1-u_mat(j-1,i+1)*v2)/(2*delta_y)) +...
%            mu*((u_mat(j,i+1)-2*u_mat(j,i)+u_mat(j,i-1))/delta_x^2 + (u_mat(j+1,i+1)-2*u_mat(j,i+1)+u_mat(j-1,i+1))/delta_y^2);
%
%        rhou(j,i+1) = rhou(j,i+1) + A*delta_t - delta_t/delta_x*(p_mat(j,i)-p_mat(j,i-1));
%
%        % Obtain velocities
%        u_mat(j,i) = rhou(j,i)/rho;
%
%    end

    % Obtain rho*v
    for i = 1:NX
        for j = 1:NY-1
            u1 = (u_mat(j,i+1)+u_mat(j+1,i+1))/2;
            u2 = (u_mat(j,i)+u_mat(j+1,i))/2;
            B = -rho*((v_mat(j+1,i+2)*u1-v_mat(j+1,i)*u2)/(2*delta_x)+(v_mat(j+2,i+1)^2-v_mat(j,i+1)^2)/(2*delta_y)) +...
                mu*((v_mat(j+1,i+2)-2*v_mat(j+1,i+1)-v_mat(j+1,i))/delta_x^2 + (v_mat(j+2,i+1)-2*v_mat(j+1,i+1)+v_mat(j,i+1))/delta_y^2);
            rhov(j+1,i+1) = rhov(j+1,i+1) + B*delta_t - delta_t/delta_y*(p_mat(j+1,i)-p_mat(j,i));

            % Obtain velocities
            v_mat(j+1,i+1) = rhov(j+1,i+1)/rho;
        end
    end

    % Apply boundary conditions
    u_mat(:,1) = u_mat(:,2); % Inflow
    u_mat(:,NXu) = u_mat(:,NXu-1); % outflow
    v_mat(:,NXv) = v_mat(:,NXv-1); % Outflow
%    v_mat(1,:) = -v_mat(2,:); % Bottom wall
%    v_mat(NY,:) = -v_mat(NYv-1,:); % Top wall

    % Keep p_mat to use later
    p_star = p_mat;

    % Solve for the pressure correction
    for loop2 = 1:iter2
        p_old = p_mat;

        % Calculate the pressure at every internal point
        for i = 2:NX-1
            for j = 2:NY-1
                a = 2*(delta_t/delta_x^2+delta_t/delta_y^2);
                b = -delta_t/delta_x^2;
                c = -delta_t/delta_y^2;
                d = (rhou(j,i+1)-rhou(j,i))/delta_x + (rhov(j+1,i+1)-rhov(j,i+1))/delta_y;
                p_mat(j,i) = -(b*p_mat(j,i+1)+b*p_mat(j,i-1)+c*p_mat(j+1,i)+c*p_mat(j-1,i)+d)/a;
            end
        end

        % Check for convergence
        if max(max(abs(p_mat-p_old))) <= tol2
            disp(['Pressure - convergence (',num2str(loop2),' iterations)'])
            break
        end

    end

    % Calculate the new values of pressure for the next iteration
    p_mat = p_star + alpha_p*p_mat;

    % Apply pressure boundary conditions
%    p_mat(1,2:end-1) = p_mat(2,2:end-1); % Bottom wall
%    p_mat(end,2:end-1) = p_mat(end-1,2:end-1); % Top wall
%    p_mat(2:end-1,1) = p_mat(2:end-1,2); % Left wall
%    p_mat(2:end-1,end) = p_mat(2:end-1,end-1); % Right wall

end

% Configure a finer grid (p)
[Xp,Yp] = meshgrid(0:delta_x:L,0:delta_y:D); % Original mesh
[Xpn,Ypn] = meshgrid(0:delta_x/2:L,0:delta_y/2:D); % Refined mesh
p_mat_n = interp2(Xp,Yp,p_mat,Xpn,Ypn); % Refined matrix

% Configure a finer grid (u)
[Xu,Yu] = meshgrid(0:delta_x:(L+delta_x),0:delta_y:D); % Original mesh
[Xun,Yun] = meshgrid(0:delta_x/2:(L+delta_x),0:delta_y/2:D); % Refined mesh
u_mat_n = interp2(Xu,Yu,u_mat,Xun,Yun); % Refined matrix
u_mat_n = u_mat_n(:,2:end-1); % Ignore points outside of the domain

% Configure a finer grid (v)
[Xv,Yv] = meshgrid(0:delta_x:(L+2*delta_x),0:delta_y:(D+delta_y)); % Original mesh
[Xvn,Yvn] = meshgrid(0:delta_x/2:(L+2*delta_x),0:delta_y/2:(D+delta_y)); % Refined mesh
v_mat_n = interp2(Xv,Yv,v_mat,Xvn,Yvn); % Refined matrix
v_mat_n = v_mat_n(2:end-1,3:end-2); % Ignore points outside of the domain

% Obtain velocity resultants
%V_mat = sqrt(u_mat_n.^2 + v_mat_n.^2);

% Make quiver plot
figure(1),clf
%quiver(Xpn,Ypn,flip(u_mat_n),flip(v_mat_n),'r','showarrowhead','off','marker','>'),grid on
quiver(Xpn,Ypn,flip(u_mat_n),flip(v_mat_n),'r'),grid on


% Obtain dimensions for a complete grid
%NX2 = 2*NX - 1;
%NY2 = 2*NY - 1;

% Plot u distribution along a vertical line in the duct
n = 10; % Number of cell along x direction
figure(2),clf
plot(flip(u_mat(:,n)),[0:delta_y:D]),grid on
xlabel('u'),ylabel('h')

% Plot v distribution in the selected point
figure(3),clf
plot([1:iter1],v_his),grid on









%% Interpolate all results into grids of equal sizes
%function sol = intp_mat(mat)
%    x = size(mat,1);
%    y = size(mat,2);
%    x2 = 2*x-1;
%    y2 = 2*y-1;
%    sol = zeros(x2,y2);
%
%    % Reorganize initial values
%    idx = 1;
%    for i = 1:2:x2
%        idy = 1;
%        for j = 1:2:y2
%            sol(i,j) = mat(idx,idy);
%            idy = idy + 1;
%        end
%        idx = idx + 1;
%    end
%    % Interpolate
%end




%    v1 = (v_mat(j+1,i)+v_mat(j+1,i+1))/2;
%    v2 = (v_mat(j,i)+v_mat(j,i+1))/2;
%    A = -rho*((u_mat(j,i+2)^2-u_mat(j,i)^2)/(2*delta_x) + (u_mat(j+1,i+1)*v1-u_mat(j-1,i+1)*v2)/(2*delta_y)) +...
%        mu*((u_mat(j,i+2)-2*u_mat(j,i+1)+u_mat(j,i))/delta_x^2 + (u_mat(j+1,i+1)-2*u_mat(j,i+1)+u_mat(j-1,i+1))/delta_y^2);
%
%    u1 = (u_mat(j,i+1)+u_mat(j+1,i+1))/2;
%    u2 = (u_mat(j,i)+u_mat(j+1,i))/2;
%    B = -rho*((v_mat(j+1,i+1)*u1-v_mat(j+1,i-1)*u2)/(2*delta_x)+(v_mat(j+2,i)^2-v_mat(j,i)^2)/(2*delta_y)) +...
%        mu*((v_mat(j+1,i+1)-2*v_mat(j+1,i)-v_mat(j+1,i-1))/delta_x^2 + (v_mat(j+2,i)-2*v_mat(j+1,i)+v_mat(j,i))/delta_y^2);
%
%    rhou(j,i+1) = rhou(j,i+1) + A*delta_t - delta_t/delta_x*(p_mat(j,i+1)-p_mat(j,i));
%    rhov(j+1,i) = rhov(j+1,i) + B*delta_t - delta_t/delta_y*(p_mat(j+1,i)-p_mat(j,i));



%%    a*p_mat(j,i)+b*p_mat(j,i+1)+b*p_mat(j,i-1)+c*p_mat(j+1,i)+c*p_mat(j-1,i)+d=0;
%    a = 2*(delta_t/delta_x^2+delta_t/delta_y^2);
%    b = -delta_t/delta_x^2;
%    c = -delta_t/delta_y^2;
%    d = (rhou(j,i+1)-rhou(j,i))/delta_x + (rhov(j+1,i)-rhov(j,i))/delta_y;
%    p_mat(j,i) = -(b*p_mat(j,i+1)+b*p_mat(j,i-1)+c*p_mat(j+1,i)+c*p_mat(j-1,i)+d)/a;
