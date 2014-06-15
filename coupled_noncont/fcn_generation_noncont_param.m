%% variables
n=5; %dim q
x_one=sym('x_one', [n 1]);
x_two=sym('x_two', [n 1]);
x=[x_one; x_two];
m=2;
%u=sym('u', [m 1]);
s=m*(2*3+1);
lam=sym('lam', [s 1]);
syms t;
taueps=sym('taueps', [4 1]);

%%parameters
a=0.730;
b=0.350;
R=0.127;
mp=21.107;  
mw=2.380;
Iw11=0.015;
Iw33=0.009;
Ip33=1.991;
ap1=0.377;
ap2=0;
gr=9.81;
T=5;

%% Q elements
Q11=mp+4*mw;
Q22=Q11;

Q44=Iw33;
Q55=Q44;
Q66=Q55;
Q77=Q66;

Q33= Ip33 + mp*(ap1^2+ap2^2) + 4*(Iw11 + mw*b^2) + 2*mw*a^2;

Q13 = -mp*(ap1*sin(x_one(3)/a) + ap2*cos(x_one(3)/a)) - 2*mw*a*sin(x_one(3)/a);
Q23 =  mp*(ap1*cos(x_one(3)/a) - ap2*sin(x_one(3)/a)) + 2*mw*a*cos(x_one(3)/a);

%% P matrix
P = [Q11   0     Q13/a   0       0;
     0     Q22   Q23/a   0       0;
     Q13/a Q23/a Q33/a^2 0       0;
     0     0     0       Q44/R^2 0;
     0     0     0       0       Q55/R^2];
 
 %% D matrix
 D = x(8)^2/a^2*([-Q23; Q13; 0; 0; 0]);
 
 %% H matrix
 
   H1=[-sin(x(3)./a) cos(x(3)./a)    0  0   0];
   H2=[-sin(x(3)./a) cos(x(3)./a)    1  0   0];
   H3=[ cos(x(3)./a) sin(x(3)./a) -b/a -1   0];
   H4=[ cos(x(3)./a) sin(x(3)./a)  b/a  0  -1];
   
  %slips
  s14=H1*x(6:10);
  s23=H2*x(6:10);
  s12=H3*x(6:10);
  s34=H4*x(6:10);
  
 %% F vector
% N forces
 N1=(mp/4+mw)*gr;
 N2=(mp/4+mw)*gr;
 N3=(mp/4+mw)*gr;
 N4=(mp/4+mw)*gr;
 
% \epsilon and \tau coefficients
 eps1=taueps(1);
 eps2=taueps(2);
 eps3=taueps(2);
 eps4=taueps(1);
 tau1=taueps(3);
 tau2=taueps(4);
 
% R coefficients
 R14=-(eps1*N1 + eps4*N4)*s14;
 R23=-(eps2*N2 + eps3*N3)*s23;
 R12=-(tau1*N1 + tau1*N2)*s12;
 R34=-(tau2*N3 + tau2*N4)*s34;
 
 H1T=[-sin(x(3)./a); cos(x(3)./a);    0;  0;   0];
   H2T=[-sin(x(3)./a); cos(x(3)./a);    1;  0;   0];
   H3T=[ cos(x(3)./a); sin(x(3)./a); -b/a; -1;   0];
   H4T=[ cos(x(3)./a); sin(x(3)./a);  b/a;  0;  -1];
 
 norm1=1;
 norm2=sqrt(2);
 norm3=sqrt(2+(b/a).^2);
 norm4=norm3;
 
 F14=R14*H1T./norm1;
 F23=R23*H2T./norm2;
 F12=R12*H3T./norm3;
 F34=R34*H4T./norm4; 
 
 F=F14+F23+F12+F34;
 
 %% Ps matrix --- control function base matrix
omega=2*pi/T;
P_vec=[1 sin(omega*t) cos(omega*t) sin(omega*2*t) cos(omega*2*t) sin(omega*3*t) cos(omega*3*t)];
P_s=blkdiag(P_vec, P_vec);

 %% control system
 B=[zeros(2, 3) eye(2)]';
g=[zeros(5,2); (1/R)*inv(P)*B];
B_lin=g;

%% output function
k=[x(1); x(2); x(3)/a];
C_lin=jacobian(k, x);
matlabFunction(k, 'file', 'sfun_k', 'vars', {x});
matlabFunction(B_lin, 'file', 'sfun_B', 'vars', {x});
matlabFunction(C_lin, 'file', 'sfun_C', 'vars', {x});
matlabFunction(g, 'file', 'sfun_g', 'vars', {x});
matlabFunction(P_s, 'file', 'sfun_P', 'vars', t);


 %% generation of A_lin and f
f=[x_two; inv(P)*(-D+F)];
A_lin=jacobian(f+ g*P_s*lam, x);
matlabFunction(f, 'file', 'sfun_f', 'vars', {x, taueps});
matlabFunction(A_lin, 'file', 'sfun_A', 'vars', {x, taueps});