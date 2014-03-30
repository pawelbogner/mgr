
emax=0.01;
kmax=200;
e=emax+1;
k=1; %iterator
r=3;
global n;
n=10;
a=0.730; %platform geometric parameter
q0=[0, 0, a*pi/2, zeros(1,2), 0, 0, 0, 0, 0]';
s=6;
%lambda=repmat([5, 1, 1]', [2 1]);
%lambda=[5 10 1 1 1 -5 10 -1 -1 -1]';
%lambda=repmat([0.5, 0.01, 0.01, 0.001, 0.001]', [2 1]);
lambda=[1, 0.1, 0.1, 0.01, 0.01, 1, 0.1, 0.1, 0.01, 0.01]';
T=20;
y_d = [10 0 pi/2 ]';

gamma=0.05;
e_tab=zeros(1, kmax);
while norm(e)>emax && k<kmax
    % tutaj trzeba pocalkowac rownanie z x i ksi
    %[t, x_ksi]=ode23(@ode_main_fun, [0, T], x_ksi);
    [t, x_ksi]=ode_run(lambda);
    %wyciagamy x i ksi
    x=x_ksi(end, 1:n);
    x=x';
    Ksi=reshape(x_ksi(end, n+1:end), n, []);
    y=sfun_k(x);
    e=y-y_d
    C=sfun_C(x);
    J=[C*Ksi];
    JP = mp_inverse(J);
    delta=JP*e;
    lambda= lambda - gamma*delta;
    [delta lambda]
    enorm=norm(e)
    e_tab(k)=enorm;
    k=k+1
%     figure(1);
%     plot(x_ksi(:,1), x_ksi(:,2));
%     axis equal
%     s=calc_slips(x_ksi(:,1:10));
%     figure(2);
%     subplot(2, 2, 1);
%     plot(t, s(:,1));
%     subplot(2, 2, 2);
%     plot(t, s(:,2));
%     subplot(2, 2, 3);
%     plot(t, s(:,3));
%     subplot(2, 2, 4);
%     plot(t, s(:,4));
%     pause;
end
figure(1);
    plot(x_ksi(:,1), x_ksi(:,2));
    axis equal
    s=calc_slips(x_ksi(:,1:10));
figure(2);
    subplot(2, 2, 1);
    plot(t, s(:,1));
    subplot(2, 2, 2);
    plot(t, s(:,2));
    subplot(2, 2, 3);
    plot(t, s(:,3));
    subplot(2, 2, 4);
    plot(t, s(:,4));
figure
    plot(e_tab);    