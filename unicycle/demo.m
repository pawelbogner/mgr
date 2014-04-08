endogen

figure(2)
plot(x_ksi(:,1), x_ksi(:,2))

figure(3)
plot(t, x_ksi(:,3)*180/pi)

figure(4)
s=calc_slips(x_ksi(:,1:8));
plot(t, s)
hold on
plot(t, 0.15*ones(1, length(t)), 'r')
plot(t, -0.15*ones(1, length(t)), 'r')

figure(5)
u=calc_u(t, lambda);
plot(t, u);
