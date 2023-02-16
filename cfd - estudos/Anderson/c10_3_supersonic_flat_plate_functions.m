1;

function delta_t = calc_delta_t(u_mat,v_mat,a_mat,rho_mat,mu_mat,delta_x,delta_y,gamma,Pr,K)
    v_prime = max(max( 4/3*mu_mat.*(gamma*mu_mat/Pr)./rho_mat ));
    delta_t_CFL = (abs(u_mat)/delta_x+abs(v_mat)/delta_y+a_mat*sqrt(delta_x^-2+delta_y^-2)+...
              2*v_prime*(delta_x^-2+delta_y^-2)).^-1;
    delta_t = min(min( K*delta_t_CFL ));
end
function mu = calc_mu(T,T1,mu1)
   mu = mu1*(T/T1).^(3/2).*(T1+110)./(T+110);
end
function k = calc_k(Pr,mu,cp)
    k = mu*cp/Pr;
end
%function [E1,E2,E3,E4] = calc_E(u,v,p,rho,tau_xx,tau_xy,Et,qx)
%    E1 = rho*u;
%    E2 = rho*u^2 + p - tau_xx;
%    E3 = rho*u*v - tau_xy;
%    E4 = (Et+p)*u - u*tau_xx - v*tau_xy + qx;
%end
%function [F1,F2,F3,F4] = calc_F(u,v,p,rho,tau_yy,tau_xy,Et,qy)
%    F1 = rho*v;
%    F2 = rho*u*v - tau_xy;
%    F3 = rho*v^2 + p - tau_yy;
%    F4 = (Et+p)*v - u*tau_xy - v*tau_yy + qy;
%end
%function tau_xx = calc_tau_xx(u1,u2,v1,v2,mu,delta_x,delta_y)
%    lambda = -2/3*mu;
%    tau_xx = lambda*((u1-u2)/delta_x+(v1-v2)/(2*delta_y)) + 2*mu*(u1-u2)/delta_x;
%end
%function tau_yy = calc_tau_yy(u1,u2,v1,v2,mu,delta_x,delta_y)
%    lambda = -2/3*mu;
%    tau_yy = lambda*((u1-u2)/(2*delta_x)+(v1-v2)/delta_y) + 2*mu*(v1-v2)/delta_y;
%end
%function tau_xy_E = calc_tau_xy_E(u1,u2,v1,v2,mu,delta_x,delta_y)
%    tau_xy_E = mu*((u1-u2)/(2*delta_y)+(v1-v2)/delta_x);
%end
%function tau_xy_F = calc_tau_xy_F(u1,u2,v1,v2,mu,delta_x,delta_y)
%    tau_xy_F = mu*((u1-u2)/delta_y+(v1-v2)/(2*delta_x));
%end
%function tau_xx = calc_tau_xx_alt(u1,u2,v1,v2,mu,delta_x,delta_y)
%    lambda = -2/3*mu;
%    tau_xx = lambda*((u1-u2)/delta_x+(v1-v2)/delta_y) + 2*mu*(u1-u2)/delta_x;
%end
%function tau_yy = calc_tau_yy_alt(u1,u2,v1,v2,mu,delta_x,delta_y)
%    lambda = -2/3*mu;
%    tau_yy = lambda*((u1-u2)/delta_x+(v1-v2)/delta_y) + 2*mu*(v1-v2)/delta_y;
%end
%function tau_xy_E = calc_tau_xy_E_alt(u1,u2,v1,v2,mu,delta_x,delta_y)
%    tau_xy_E = mu*((u1-u2)/delta_y+(v1-v2)/delta_x);
%end
%function tau_xy_F = calc_tau_xx_alt(u1,u2,v1,v2,mu,delta_x,delta_y)
%    tau_xx = mu*((u1-u2)/delta_y+(v1-v2)/delta_x);
%end
function [qx,qy] = calc_q(k,T1x,T2x,T1y,T2y,delta_x,delta_y)
    qx = -k*(T1x-T2x)/delta_x;
    qy = -k*(T1y-T2y)/delta_y;
end
function Et = calc_Et(u,v,rho,e)
    V = sqrt(u.^2+v.^2);
    Et = rho.*(e+V.^2/2);
end
function [u,v,rho,p,T,e] = air_from_U(U1,U2,U3,U4,cv,R)
    rho = U1;
    u = U2./U1;
    v = U3./U1;
    Et = U4;

%    e = U4./U1-(u.^2+v.^2)/2;
%    T = e/cv;
%    p = rho.*R.*T;

    V = sqrt(u.^2+v.^2);
    e = Et./rho-V.^2/2;
    T = e/cv;
    p = rho.*R.*T;
end
function [U1,U2,U3,U4] = U_from_air(u,v,rho,e)
    U1 = rho;
    U2 = rho.*u;
    U3 = rho.*v;
    V = sqrt(u.^2+v.^2);
    U4 = rho.*(e+V.^2/2);
end
function [u,v,T,p,rho,e,mu] = boundaries(u,v,T,p,rho,e,mu,u1,v1,T1,p1,rho1,mu1,cv,R,Tw,NX,NY)
    % Case 1 - Leading edge (range: point 1,1 only)
    u(1,1) = 0;
    v(1,1) = 0;
    T(1,1) = T1;
    p(1,1) = p1;
    rho(1,1) = rho1;
    e(1,1) = cv*T1;
    mu(1,1) = mu1;
    % Case 2 - Right and top boundaries (range: (2:NY,1) and (NY,2:NX))
    u(2:end,1) = repmat(u1,NY-1,1); % Left boundary
    u(end,2:end) = repmat(u1,1,NX-1); % Top boundary
    v(2:end,1) = zeros(NY-1,1); % Left boundary
    v(end,2:end) = zeros(1,NX-1); % Top boundary
    T(2:end,1) = repmat(T1,NY-1,1); % Left boundary
    T(end,2:end) = repmat(T1,1,NX-1); % Top boundary
    p(2:end,1) = repmat(p1,NY-1,1); % Left boundary
    p(end,2:end) = repmat(p1,1,NX-1); % Top boundary
    rho(2:end,1) = repmat(rho1,NY-1,1); % Left boundary
    rho(end,2:end) = repmat(rho1,1,NX-1); % Top boundary
    e(2:end,1) = repmat(cv*T1,NY-1,1); % Left boundary
    e(end,2:end) = repmat(cv*T1,1,NX-1); % Top boundary
    mu(2:end,1) = repmat(mu1,NY-1,1); % Left boundary
    mu(end,2:end) = repmat(mu1,1,NX-1); % Top boundary
    % Case 4 - Outlet (Left) (range: (2:NY-1,NX))
    % (case 4 affects the last point of the plate, so it comes before case 3)
    % (this could be changed if the boundary region ranges were changed)
    u(2:NY-1,NX) = 2*u(2:NY-1,NX-1) - u(2:NY-1,NX-2);
    v(2:NY-1,NX) = 2*v(2:NY-1,NX-1) - v(2:NY-1,NX-2);
    T(2:NY-1,NX) = 2*T(2:NY-1,NX-1) - T(2:NY-1,NX-2);
    p(2:NY-1,NX) = 2*p(2:NY-1,NX-1) - p(2:NY-1,NX-2);
    rho(2:NY-1,NX) = p(2:NY-1,NX)./(R*T(2:NY-1,NX));
    e(2:NY-1,NX) = cv*T(2:NY-1,NX);
    mu(2:NY-1,NX) = calc_mu(T(2:NY-1,NX),T1,mu1);
    % Case 3 - Wall (bottom) (range: (1,2:NX))
    u(1,2:end) = zeros(1,NX-1);
    v(1,2:end) = zeros(1,NX-1);
    T(1,2:end) = repmat(Tw,1,NX-1);
    p(1,2:end) = 2*p(2,2:end) - p(3,2:end);
    rho(1,2:end) = p(1,2:end)./(R*T(1,2:end));
%    rho(1,2:end) = zeros(1,NX-1);
    e(1,2:end) = repmat(cv*Tw,1,NX-1);
    mu(1,2:end) = repmat(calc_mu(Tw,T1,mu1),1,NX-1);
end

%function [E1,E2,E3,E4,F1,F2,F3,F4] = E_F_borders(u_mat,v_mat,p_mat,rho_mat,T_mat,e_mat,mu_mat,cp,Pr,delta_x,delta_y,NX,NY)
%    j = 1;
%    for i = 1:NX
%        if i < NX
%            % Calculate tau terms
%            tau_xx = calc_tau_xx_alt(u_mat(j,i+1),u_mat(j,i),v_mat(j+1,i),v_mat(j,i),mu_mat(j,i),delta_x,delta_y);
%            tau_yy = calc_tau_yy_alt(u_mat(j,i+1),u_mat(j,i),v_mat(j+1,i),v_mat(j,i),mu_mat(j,i),delta_x,delta_y);
%            tau_xy_E = calc_tau_xy_E_alt(u_mat(j+1,i),u_mat(j,i),v_mat(j,i+1),v_mat(j,i),mu_mat(j,i),delta_x,delta_y);
%            tau_xy_F = tau_xy_E;
%
%            % Calculate qx and qy
%            k = calc_k(Pr,mu_mat(j,i),cp);
%            [qx,qy] = calc_q(k,T_mat(j,i+1),T_mat(j,i),T_mat(j+1,i),T_mat(j,i),delta_x,delta_y);
%        else
%            % Calculate tau terms
%            tau_xx = calc_tau_xx_alt(u_mat(j,i),u_mat(j,i-1),v_mat(j+1,i),v_mat(j,i),mu_mat(j,i),delta_x,delta_y);
%            tau_yy = calc_tau_yy_alt(u_mat(j,i),u_mat(j,i-1),v_mat(j+1,i),v_mat(j,i),mu_mat(j,i),delta_x,delta_y);
%            tau_xy_E = calc_tau_xy_E_alt(u_mat(j+1,i),u_mat(j,i),v_mat(j,i),v_mat(j,i-1),mu_mat(j,i),delta_x,delta_y);
%            tau_xy_F = tau_xy_E;
%
%            % Calculate qx and qy
%            k = calc_k(Pr,mu_mat(j,i),cp);
%            [qx,qy] = calc_q(k,T_mat(j,i),T_mat(j,i-1),T_mat(j+1,i),T_mat(j,i),delta_x,delta_y);
%        end
%        % Lower boundary
%        Et = calc_Et(u_mat(j,i),v_mat(j,i),rho_mat(j,i),e_mat(j,i));
%        [E1(1,i),E2(1,i),E3(1,i),E4(1,i)] = calc_E(u_mat(1,i),v_mat(1,i),p_mat(1,i),rho_mat(1,i),tau_xx,tau_xy_E,Et,qx);
%        [F1(1,i),F2(1,i),F3(1,i),F4(1,i)] = calc_F(u_mat(1,i),v_mat(1,i),p_mat(1,i),rho_mat(1,i),tau_yy,tau_xy_F,Et,qy);
%    end
%    % Obtain values of E and F along the border (upper)
%    j = NY;
%    for i = 1:NX
%        if i < NX
%            % Calculate tau terms
%            tau_xx = calc_tau_xx_alt(u_mat(j,i+1),u_mat(j,i),v_mat(j,i),v_mat(j-1,i),mu_mat(j,i),delta_x,delta_y);
%            tau_yy = calc_tau_yy_alt(u_mat(j,i+1),u_mat(j,i),v_mat(j,i),v_mat(j-1,i),mu_mat(j,i),delta_x,delta_y);
%            tau_xy_E = calc_tau_xy_E_alt(u_mat(j,i),u_mat(j-1,i),v_mat(j,i+1),v_mat(j,i),mu_mat(j,i),delta_x,delta_y);
%            tau_xy_F = tau_xy_E;
%
%            % Calculate qx and qy
%            k = calc_k(Pr,mu_mat(j,i),cp);
%            [qx,qy] = calc_q(k,T_mat(j,i+1),T_mat(j,i),T_mat(j,i),T_mat(j-1,i),delta_x,delta_y);
%        else
%            % Calculate tau terms
%            tau_xx = calc_tau_xx_alt(u_mat(j,i),u_mat(j,i-1),v_mat(j,i),v_mat(j-1,i),mu_mat(j,i),delta_x,delta_y);
%            tau_yy = calc_tau_yy_alt(u_mat(j,i),u_mat(j,i-1),v_mat(j,i),v_mat(j-1,i),mu_mat(j,i),delta_x,delta_y);
%            tau_xy_E = calc_tau_xy_E_alt(u_mat(j,i),u_mat(j-1,i),v_mat(j,i),v_mat(j,i-1),mu_mat(j,i),delta_x,delta_y);
%            tau_xy_F = tau_xy_E;
%
%            % Calculate qx and qy
%            k = calc_k(Pr,mu_mat(j,i),cp);
%            [qx,qy] = calc_q(k,T_mat(j,i),T_mat(j,i-1),T_mat(j,i),T_mat(j-1,i),delta_x,delta_y);
%        end
%        % Upper boundary
%        Et = calc_Et(u_mat(j,i),v_mat(j,i),rho_mat(j,i),e_mat(j,i));
%        [E1(j,i),E2(j,i),E3(j,i),E4(j,i)] = calc_E(u_mat(j,i),v_mat(j,i),p_mat(j,i),rho_mat(j,i),tau_xx,tau_xy_E,Et,qx);
%        [F1(j,i),F2(j,i),F3(j,i),F4(j,i)] = calc_F(u_mat(j,i),v_mat(j,i),p_mat(j,i),rho_mat(j,i),tau_yy,tau_xy_F,Et,qy);
%    end
%    % Obtain values of E and F along the border (left)
%    i = 1;
%    for j = 2:NY-1
%        % Calculate tau terms
%        tau_xx = calc_tau_xx_alt(u_mat(j,i+1),u_mat(j,i),v_mat(j+1,i),v_mat(j,i),mu_mat(j,i),delta_x,delta_y);
%        tau_yy = calc_tau_yy_alt(u_mat(j,i+1),u_mat(j,i),v_mat(j+1,i),v_mat(j,i),mu_mat(j,i),delta_x,delta_y);
%        tau_xy_E = calc_tau_xy_E_alt(u_mat(j+1,i),u_mat(j,i),v_mat(j,i+1),v_mat(j,i),mu_mat(j,i),delta_x,delta_y);
%        tau_xy_F = tau_xy_E;
%
%        % Calculate qx and qy
%        k = calc_k(Pr,mu_mat(j,i),cp);
%        [qx,qy] = calc_q(k,T_mat(j,i+1),T_mat(j,i),T_mat(j+1,i),T_mat(j,i),delta_x,delta_y);
%
%        % Left boundary
%        Et = calc_Et(u_mat(j,i),v_mat(j,i),rho_mat(j,i),e_mat(j,i));
%        [E1(j,i),E2(j,i),E3(j,i),E4(j,i)] = calc_E(u_mat(j,i),v_mat(j,i),p_mat(j,i),rho_mat(j,i),tau_xx,tau_xy_E,Et,qx);
%        [F1(j,i),F2(j,i),F3(j,i),F4(j,i)] = calc_F(u_mat(j,i),v_mat(j,i),p_mat(j,i),rho_mat(j,i),tau_yy,tau_xy_F,Et,qy);
%    end
%    % Obtain values of E and F along the border (right)
%    i = NX;
%    for j = 2:NY-1
%        % Calculate tau terms
%        tau_xx = calc_tau_xx_alt(u_mat(j,i),u_mat(j,i-1),v_mat(j+1,i),v_mat(j,i),mu_mat(j,i),delta_x,delta_y);
%        tau_yy = calc_tau_yy_alt(u_mat(j,i),u_mat(j,i-1),v_mat(j+1,i),v_mat(j,i),mu_mat(j,i),delta_x,delta_y);
%        tau_xy_E = calc_tau_xy_E_alt(u_mat(j+1,i),u_mat(j,i),v_mat(j,i),v_mat(j,i-1),mu_mat(j,i),delta_x,delta_y);
%        tau_xy_F = tau_xy_E;
%
%        % Calculate qx and qy
%        k = calc_k(Pr,mu_mat(j,i),cp);
%        [qx,qy] = calc_q(k,T_mat(j,i),T_mat(j,i-1),T_mat(j+1,i),T_mat(j,i),delta_x,delta_y);
%
%        % Left boundary
%        Et = calc_Et(u_mat(j,i),v_mat(j,i),rho_mat(j,i),e_mat(j,i));
%        [E1(j,i),E2(j,i),E3(j,i),E4(j,i)] = calc_E(u_mat(j,i),v_mat(j,i),p_mat(j,i),rho_mat(j,i),tau_xx,tau_xy_E,Et,qx);
%        [F1(j,i),F2(j,i),F3(j,i),F4(j,i)] = calc_F(u_mat(j,i),v_mat(j,i),p_mat(j,i),rho_mat(j,i),tau_yy,tau_xy_F,Et,qy);
%    end
%
%
%end

