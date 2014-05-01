function u = calc_u(t, lambda)
N=size(t, 1);
u1=zeros(1,N);
u2=zeros(1,N);
u=[u1' u2'];
for i=1:N
    u(i,:)=(sfun_P(t(i))*lambda)';
end