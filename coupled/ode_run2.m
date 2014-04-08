function [tout, xout] = ode_run2(lambdas)

tstart=0;
tfinal=4;
x0=[0; 0; 0.730*pi/2; zeros(7,1)];

odefcn=@(t,x) sfun_f(x)+sfun_g(x)*sfun_P(t)*lambdas;

[tout, xout] = ode45(odefcn, [tstart, tfinal], x0);