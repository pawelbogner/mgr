%% parameters
m=2;
R=0.3;
Iph=m*R^2/4;
Ith=m*R^2/2;

gr=9.81;

T=4;

epss=1;
taus=10;
epsd=0.5;
taud=0.3;

n=4;
m=2;

q=sym('q', [n 1]);
qd=sym('qd', [n 1]);
x=[q; qd];
syms t;
s=m*(2*1+1);
lam=sym('lam', [s 1]);

%% inertia matrix
M=diag([m, m, Iph, Ith]);

%% control input matrix
B=[zeros(2,2); eye(2)];

%% pfaffian matrix
H=[-sin(q(3))  cos(q(3)) 0  0;
    cos(q(3))  sin(q(3)) 0 -R];

%% slips
s1=H(1,:)*qd;
s2=H(2,:)*qd;

%% friction reaction forces
R1s=-epss*m*gr*s1;
R2s=-taus*m*gr*s2;
R1d=-epsd*m*gr*s1;
R2d=-taud*m*gr*s2;

%% friction
F0=R1s*H(1,:)'/norm(H(1,:))+R2s*H(2,:)'/norm(H(2,:));
F1=R1d*H(1,:)'/norm(H(1,:))+R2s*H(2,:)'/norm(H(2,:));
F2=R1s*H(1,:)'/norm(H(1,:))+R2d*H(2,:)'/norm(H(2,:));
F3=R1d*H(1,:)'/norm(H(1,:))+R2d*H(2,:)'/norm(H(2,:));

%% Ps matrix --- control function base matrix
omega=2*pi/T;
P_vec=[1 sin(omega*t) cos(omega*t)];
P_s=blkdiag(P_vec, P_vec);

%% output function
k=[q(1); q(2); q(3)];

%% system equation
f0=[qd; M\F0];
f1=[qd; M\F1];
f2=[qd; M\F2];
f3=[qd; M\F3];
g=[zeros(4, 2); M\B];

A_lin0=jacobian(f0+ g*P_s*lam, x);
A_lin1=jacobian(f1+ g*P_s*lam, x);
A_lin2=jacobian(f2+ g*P_s*lam, x);
A_lin3=jacobian(f3+ g*P_s*lam, x);
B_lin=g;
C_lin=jacobian(k, x);

%% generating m-files
matlabFunction(k, 'file', 'sfun_k', 'vars', {x});
matlabFunction(f0, 'file', 'sfun_f0', 'vars', {x});
matlabFunction(f1, 'file', 'sfun_f1', 'vars', {x});
matlabFunction(f2, 'file', 'sfun_f2', 'vars', {x});
matlabFunction(f3, 'file', 'sfun_f3', 'vars', {x});
matlabFunction(g, 'file', 'sfun_g');
matlabFunction(P_s, 'file', 'sfun_P', 'vars', t);

matlabFunction(A_lin0, 'file', 'sfun_A0', 'vars', {x, lam});
matlabFunction(A_lin1, 'file', 'sfun_A1', 'vars', {x, lam});
matlabFunction(A_lin2, 'file', 'sfun_A2', 'vars', {x, lam});
matlabFunction(A_lin3, 'file', 'sfun_A3', 'vars', {x, lam});
matlabFunction(B_lin, 'file', 'sfun_B');
matlabFunction(C_lin, 'file', 'sfun_C', 'vars', {x});