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
function Et = calc_Et(u,v,rho,e)
    V = sqrt(u.^2+v.^2);
    Et = rho.*(e+V.^2/2);
end
function [u,v,rho,p,T,e] = air_from_U(U1,U2,U3,U4,cv,R)
    rho = U1;
    u = U2./U1;
    v = U3./U1;
    Et = U4;
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
function [E1,E2,E3,E4] = E_from_air(u,v,p,rho,e,tau_xx,tau_xy,qx)
    E1 = rho.*u;
    E2 = rho.*u.^2 + p - tau_xx;
    E3 = rho.*u.*v - tau_xy;
    Et = calc_Et(u,v,rho,e);
    E4 = (Et+p).*u - u.*tau_xx - v.*tau_xy + qx;
end
function [F1,F2,F3,F4] = F_from_air(u,v,p,rho,e,tau_yy,tau_xy,qy)
    F1 = rho.*v;
    F2 = rho.*u.*v - tau_xy;
    F3 = rho.*v.^2 + p - tau_yy;
    Et = calc_Et(u,v,rho,e);
    F4 = (Et+p).*v - u.*tau_xy - v.*tau_yy + qy;
end
function [u,v,T,p,rho,e,mu] = boundaries(u,v,T,p,rho,e,mu,u1,v1,T1,p1,rho1,mu1,cv,R,Tw,NX,NY,LEx)
    % Case 1 - Leading edge (range: points (1,LEx) and (NY,LEx) only)
    u(1,LEx) = 0;
    v(1,LEx) = 0;
    T(1,LEx) = T1;
    p(1,LEx) = p1;
    rho(1,LEx) = rho1;
    e(1,LEx) = cv*T1;
    mu(1,LEx) = mu1;
    u(NY,LEx) = 0;
    v(NY,LEx) = 0;
    T(NY,LEx) = T1;
    p(NY,LEx) = p1;
    rho(NY,LEx) = rho1;
    e(NY,LEx) = cv*T1;
    mu(NY,LEx) = mu1;
    % Case 2 - Right boundary (range: (1:NY,1))
    u(1:end,1) = repmat(u1,NY,1);
    v(1:end,1) = zeros(NY,1);
    T(1:end,1) = repmat(T1,NY,1);
    p(1:end,1) = repmat(p1,NY,1);
    rho(1:end,1) = repmat(rho1,NY,1);
    e(1:end,1) = repmat(cv*T1,NY,1);
    mu(1:end,1) = repmat(mu1,NY,1);
    % Case 2 - Bottom boundary (range: (1,2:LEx))
    u(1,2:LEx-1) = repmat(u1,1,LEx-2);
    v(1,2:LEx-1) = zeros(1,LEx-2);
    T(1,2:LEx-1) = repmat(T1,1,LEx-2);
    p(1,2:LEx-1) = repmat(p1,1,LEx-2);
    rho(1,2:LEx-1) = repmat(rho1,1,LEx-2);
    e(1,2:LEx-1) = repmat(cv*T1,1,LEx-2);
    mu(1,2:LEx-1) = repmat(mu1,1,LEx-2);
    % Case 2 - Top boundary (range: (NY,2:LEx))
    u(end,2:LEx-1) = repmat(u1,1,LEx-2);
    v(end,2:LEx-1) = zeros(1,LEx-2);
    T(end,2:LEx-1) = repmat(T1,1,LEx-2);
    p(end,2:LEx-1) = repmat(p1,1,LEx-2);
    rho(end,2:LEx-1) = repmat(rho1,1,LEx-2);
    e(end,2:LEx-1) = repmat(cv*T1,1,LEx-2);
    mu(end,2:LEx-1) = repmat(mu1,1,LEx-2);
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
    % Case 3 - Wall (bottom) (range: (1,LEx+1:NX))
    u(1,LEx+1:end) = zeros(1,NX-LEx);
    v(1,LEx+1:end) = zeros(1,NX-LEx);
    T(1,LEx+1:end) = repmat(Tw,1,NX-LEx);
    p(1,LEx+1:end) = 2*p(2,LEx+1:end) - p(3,LEx+1:end);
    rho(1,LEx+1:end) = p(1,LEx+1:end)./(R*T(1,LEx+1:end));
    e(1,LEx+1:end) = repmat(cv*Tw,1,NX-LEx);
    mu(1,LEx+1:end) = repmat(calc_mu(Tw,T1,mu1),1,NX-LEx);
    % Case 3 - Wall (top) (range: (NY,LEx+1:NX))
    u(end,LEx+1:end) = zeros(1,NX-LEx);
    v(end,LEx+1:end) = zeros(1,NX-LEx);
    T(end,LEx+1:end) = repmat(Tw,1,NX-LEx);
    p(end,LEx+1:end) = 2*p(end-1,LEx+1:end) - p(end-2,LEx+1:end);
    rho(end,LEx+1:end) = p(1,LEx+1:end)./(R*T(1,LEx+1:end));
    e(end,LEx+1:end) = repmat(cv*Tw,1,NX-LEx);
    mu(end,LEx+1:end) = repmat(calc_mu(Tw,T1,mu1),1,NX-LEx);
end
function tau_xx = calc_tau_xx_P(u,v,mu,delta_x,delta_y,NX,NY)
    % Predictor step only
    % x derivatives: backwards differences
    % y derivatives: central differences
    lambda = -2/3*mu;
    tau_xx = zeros(NY,NX);
    for i = 2:NX
        for j = 2:NY-1
            if (i==1&&j==1) || (i==1&&j==NY) || (i==NX&&j==1) || (i==NX&&j==NY)
                % Ignore corners of the domain
                continue
            end
            tau_xx(j,i) = lambda(j,i)*((u(j,i)-u(j,i-1))/delta_x+(v(j+1,i)-v(j-1,i))/(2*delta_y)) + 2*mu(j,i)*(u(j,i)-u(j,i-1))/delta_x;
        end
    end
end
function tau_yy = calc_tau_yy_P(u,v,mu,delta_x,delta_y,NX,NY)
    % Predictor step only
    % x derivatives: central differences
    % y derivatives: backwards differences
    lambda = -2/3*mu;
    tau_yy = zeros(NY,NX);
    for i = 2:NX-1
            for j = 2:NY
            if (i==1&&j==1) || (i==1&&j==NY) || (i==NX&&j==1) || (i==NX&&j==NY)
                % Ignore corners of the domain
                continue
            end
            tau_yy(j,i) = lambda(j,i)*((u(j,i+1)-u(j,i-1))/(2*delta_x)+(v(j,i)-v(j-1,i))/delta_y) + 2*mu(j,i)*(v(j,i)-v(j-1,i))/delta_y;
        end
    end
end
function tau_xy_E = calc_tau_xy_EP(u,v,mu,delta_x,delta_y,NX,NY)
    % Predictor step only, for E terms
    % x derivatives: backwards differences
    % y derivatives: central differences
    tau_xy_E = zeros(NY,NX);
    for i = 2:NX
        for j = 2:NY-1
            if (i==1&&j==1) || (i==1&&j==NY) || (i==NX&&j==1) || (i==NX&&j==NY)
                % Ignore corners of the domain
                continue
            end
            tau_xy_E(j,i) = mu(j,i)*((u(j+1,i)-u(j-1,i))/(2*delta_y)+(v(j,i)-v(j,i-1))/delta_x);
        end
    end
end
function tau_xy_F = calc_tau_xy_FP(u,v,mu,delta_x,delta_y,NX,NY)
    % Predictor step only, for F terms
    % x derivatives: central differences
    % y derivatives: backwards differences
    tau_xy_F = zeros(NY,NX);
    for i = 2:NX-1
        for j = 2:NY
            if (i==1&&j==1) || (i==1&&j==NY) || (i==NX&&j==1) || (i==NX&&j==NY)
                % Ignore corners of the domain
                continue
            end
            tau_xy_F(j,i) = mu(j,i)*((u(j,i)-u(j-1,i))/delta_y+(v(j,i+1)-v(j,i-1))/(2*delta_x));
        end
    end
end
function qx = calc_qx_P(T,k,delta_x,NX,NY)
    % Predictor step only
    % x derivative: backward difference
    qx = zeros(NY,NX);
    for i = 2:NX
        for j = 2:NY-1
            qx(j,i) = -k(j,i)*(T(j,i)-T(j,i-1))/delta_x;
        end
    end
end
function qy = calc_qy_P(T,k,delta_y,NX,NY)
    % Predictor step only
    % y derivative: backward difference
    qy = zeros(NY,NX);
    for i = 2:NX-1
        for j = 2:NY
            qy(j,i) = -k(j,i)*(T(j,i)-T(j-1,i))/delta_y;
        end
    end
end
function tau_xx = calc_tau_xx_C(u,v,mu,delta_x,delta_y,NX,NY)
    % Corrector step only
    % x derivatives: backwards differences
    % y derivatives: central differences
    lambda = -2/3*mu;
    tau_xx = zeros(NY,NX);
    for i = 1:NX-1
        for j = 2:NY-1
            if (i==1&&j==1) || (i==1&&j==NY) || (i==NX&&j==1) || (i==NX&&j==NY)
                % Ignore corners of the domain
                continue
            end
            tau_xx(j,i) = lambda(j,i)*((u(j,i+1)-u(j,i))/delta_x+(v(j+1,i)-v(j-1,i))/(2*delta_y)) + 2*mu(j,i)*(u(j,i+1)-u(j,i))/delta_x;
        end
    end
end
function tau_yy = calc_tau_yy_C(u,v,mu,delta_x,delta_y,NX,NY)
    % Corrector step only
    % x derivatives: central differences
    % y derivatives: backwards differences
    lambda = -2/3*mu;
    tau_yy = zeros(NY,NX);
    for i = 2:NX-1
        for j = 1:NY-1
            if (i==1&&j==1) || (i==1&&j==NY) || (i==NX&&j==1) || (i==NX&&j==NY)
                % Ignore corners of the domain
                continue
            end
            tau_yy(j,i) = lambda(j,i)*((u(j,i+1)-u(j,i-1))/(2*delta_x)+(v(j+1,i)-v(j,i))/delta_y) + 2*mu(j,i)*(v(j+1,i)-v(j,i))/delta_y;
        end
    end
end
function tau_xy_E = calc_tau_xy_EC(u,v,mu,delta_x,delta_y,NX,NY)
    % Corrector step only, for E terms
    % x derivatives: backwards differences
    % y derivatives: central differences
    tau_xy_E = zeros(NY,NX);
    for i = 1:NX-1
        for j = 2:NY-1
            if (i==1&&j==1) || (i==1&&j==NY) || (i==NX&&j==1) || (i==NX&&j==NY)
                % Ignore corners of the domain
                continue
            end
            tau_xy_E(j,i) = mu(j,i)*((u(j+1,i)-u(j-1,i))/(2*delta_y)+(v(j,i+1)-v(j,i))/delta_x);
        end
    end
end
function tau_xy_F = calc_tau_xy_FC(u,v,mu,delta_x,delta_y,NX,NY)
    % Corrector step only, for F terms
    % x derivatives: central differences
    % y derivatives: backwards differences
    tau_xy_F = zeros(NY,NX);
    for i = 2:NX-1
        for j = 1:NY-1
            if (i==1&&j==1) || (i==1&&j==NY) || (i==NX&&j==1) || (i==NX&&j==NY)
                % Ignore corners of the domain
                continue
            end
            tau_xy_F(j,i) = mu(j,i)*((u(j+1,i)-u(j,i))/delta_y+(v(j,i+1)-v(j,i-1))/(2*delta_x));
        end
    end
end
function qx = calc_qx_C(T,k,delta_x,NX,NY)
    % Corrector step only
    % x derivative: forward difference
    qx = zeros(NY,NX);
    for i = 1:NX-1
        for j = 2:NY-1
            qx(j,i) = -k(j,i)*(T(j,i+1)-T(j,i))/delta_x;
        end
    end
end
function qy = calc_qy_C(T,k,delta_y,NX,NY)
    % Predictor step only
    % y derivative: forward difference
    qy = zeros(NY,NX);
    for i = 2:NX-1
        for j = 1:NY-1
            qy(j,i) = -k(j,i)*(T(j+1,i)-T(j,i))/delta_y;
        end
    end
end
function status = convergence(rho,rho_old,tol)
    val = abs((rho-rho_old)./rho);
    if max(max(val)) <= tol
        status = 1;
    else
        status = 0;
    end
end
function plot_M_contours(u_mat,v_mat,a_mat)
    figure(1),clf
    M_mat = sqrt(u_mat.^2+v_mat.^2)./a_mat;
    contourf(M_mat)
    title('Mach number')
    colorbar
    colormap turbo
end
