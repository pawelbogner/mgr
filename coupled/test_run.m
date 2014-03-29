la=[0.05 0.1 0 0 0 0.1 0.1 0 0 0]';
[t, x]=ode_run(la);
figure(1);
plot(x(:,1), x(:,2));
axis equal
figure(2);
s=calc_slips(x);
subplot(2, 2, 1);
plot(t, s(:,1));
subplot(2, 2, 2);
plot(t, s(:,2));
subplot(2, 2, 3);
plot(t, s(:,3));
subplot(2, 2, 4);
plot(t, s(:,4));
figure(3);
plot(t, x(:,3)/pi*180/a);