function [tout, xout] = ode_run(lambdas)

tstart=0;
tfinal=10;
%x0=[0; 0; 0.730*pi/2; zeros(7,1)]; %% dim Ksi= n*m*s
x0=[0; 0; 0.730*pi/2; zeros(7,1); zeros(10*2*7,1)];
opt=odeset('MaxStep', 1e-2);
% epstau=[1 1 1 1 3 3 3 3]';
    % ode_fun = [f+g Ps lam; A ksi+ BP]

    odefcn=@(t,x) [sfun_f(x(1:10))+sfun_g(x(1:10))*sfun_P(t)*lambdas; reshape(sfun_A(x(1:10))*reshape(x(11:end), 10, [])+sfun_B(x(1:10))*sfun_P(t), [], 1)];
%odefcn=@(t,x) [sfun_f(x(1:10))+sfun_g(x(1:10))*sfun_P(t)*lambdas];
    
    [tout, xout] = ode15s(odefcn, [tstart, tfinal], x0, opt);