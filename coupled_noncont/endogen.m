
emax=0.01;
kmax=100;
e=emax+1;
k=1; %iterator
r=3;
global n;
n=10;
a=0.730; %platform geometric parameter
%lambda=repmat([1, 0, 0]', [2 1]);
lambda=repmat([0.5, 0.01, 0.01, 0.001, 0.001, 0.0001, 0.0001]', [2 1]);
T=4;
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
    lambda= lambda - gamma*JP*e;
    e_tab(k)=norm(e);
    k=k+1    
end

plot(e_tab);    