function [tout, xout] = start_sim(lambda)

tstart=0;
tfinal=4;
x0=[0; 0; 0.730*pi/2; zeros(7,1)];

odefcn=@(t,x) sfun_f(x)+sfun_g(x)*sfun_P(t)*lambda;

[tout, xout] = ode45(odefcn, [tstart, tfinal], x0);