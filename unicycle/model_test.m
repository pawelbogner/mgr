lambda=[0; zeros(8,1); 1; zeros(8,1)];

odefun = @(t,x) sfun_f(x)+sfun_g(x)*sfun_P(t)*lambda;

[time, traj]=ode45(odefun, [0, 4], zeros(8,1));
plot(time, traj(:,2))