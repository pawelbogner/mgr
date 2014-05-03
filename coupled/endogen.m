
emax=0.01;
kmax=200;
e=emax+1;
k=1; %iterator
r=4;
q=[pi/4 pi/5 pi/6 -pi/3]';
global n;
n=10;
a=0.730; %platform geometric parameter
%lambda=repmat([5, 1, 1]', [2 1]);
%lambda=[5 10 1 1 1 -5 10 -1 -1 -1]';
lambda=repmat([2, 0, 1, 0.5]', [2 1]);
s=8;
%lambda=[1, 0.1, 0.1, 0.01, 0.01, 1, 0.1, 0.1, 0.01, 0.01]';
T=4;
y_d = [10 0 pi/2 0.2 pi/2 pi/4 -pi/2]';

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
    y=sfun_k(x, q);
    e=y-y_d
    C=sfun_C(x);
    D=sfun_D(q);
    J=[C*Ksi D];
    %J=[C*Ksi];
    lambda_q = [lambda; q ];
    JP = mp_inverse(J);
    delta=JP*e;
    lambda_q= lambda_q - gamma*delta;
    lambda=lambda_q(1:s);
    q=lambda_q(s+1:end);
    [delta lambda_q]
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
figure(2);
    s=calc_slips(x_ksi(:,1:10));
    subplot(2, 2, 1);
    plot(t, s(:,1));
    subplot(2, 2, 2);
    plot(t, s(:,2));
    subplot(2, 2, 3);
    plot(t, s(:,3));
    subplot(2, 2, 4);
    plot(t, s(:,4));
figure(3);
    plot(t, x_ksi(:,3)*180/pi/a);
figure(4);
    plot(e_tab);    