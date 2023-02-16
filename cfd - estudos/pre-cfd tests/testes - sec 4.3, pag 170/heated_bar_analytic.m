function ans = heated_bar_analytic(x,L,T_Im,alpha,t,sum_max)
    
    if nargin == 5
        sum_max = 100;
    end
    
    sum_val = 0;
    for n = 1:sum_max
        sum_val = (-1)^n*2*T_Im/(n*pi)*sin(n*pi*x/L)*exp(-alpha*t*(n*pi/L)^2) + sum_val;
    end
    
    ans = x/L*T_Im + sum_val;

end