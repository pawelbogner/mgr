function [tout, xout] = ode_run2(lambdas)

tstart=0;
tfinal=4;
%x0=[0; 0; 0.730*pi/2; zeros(7,1)]; %% dim Ksi= n*m*s
x0=[0; 0; 0.730*pi/2; zeros(7,1)];
opt=odeset('MaxStep', 1e-2);
%epstau=[1 1 1 1 3 3 3 3]';
    % ode_fun = [f+g Ps lam; A ksi+ BP]

    odefcn=@(t,x) sfun_f(x(1:10))+sfun_g(x(1:10))*sfun_P(t)*lambdas;
%odefcn=@(t,x) [sfun_f(x(1:10))+sfun_g(x(1:10))*sfun_P(t)*lambdas];
    
    [tout, xout] = ode45(odefcn, [tstart, tfinal], x0, opt);