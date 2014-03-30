%% variables
n=5; %dim q
x_one=sym('x_one', [n 1]);
x_two=sym('x_two', [n 1]);
x=[x_one; x_two];
m=2;
%u=sym('u', [m 1]);
s=m*(2*2+1);

lam=sym('lam', [s 1]);
syms t;

%%parameters
a=0.730;
b=0.530;
R=0.127;
mp=21.107;  
mw=2.380;
Iw11=0.015;
Iw33=0.009;
Ip33=1.991;
ap1=0.377;
ap2=0;
gr=9.81;
T=20;

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
 H = [-sin(x_one(3)/a) cos(x_one(3)/a)    0  0   0;
      -sin(x_one(3)/a) cos(x_one(3)/a)    1  0   0;
       cos(x_one(3)/a) sin(x_one(3)/a) -b/a -1   0;
       cos(x_one(3)/a) sin(x_one(3)/a)  b/a  0  -1];
 
  %slips
  s14=H(1,:)*x_two;
  s23=H(2,:)*x_two;
  s12=H(3,:)*x_two;
  s34=H(4,:)*x_two;
  
 %% F vector
% N forces
 N1=(mp/4+mw)*gr;
 N2=(mp/4+mw)*gr;
 N3=(mp/4+mw)*gr;
 N4=(mp/4+mw)*gr;
 
% \epsilon and \tau coefficients
 epstau=sym('epstau', [8 1]);
 
% R coefficients
 R14=-(epstau(1)*N1 + epstau(4)*N4)*s14;
 R23=-(epstau(2)*N2 + epstau(3)*N3)*s23;
 R12=-(epstau(5)*N1 + epstau(6)*N2)*s12;
 R34=-(epstau(7)*N3 + epstau(8)*N4)*s34;
 
 F14=R14*H(1,:)'/norm(H(1,:));
 F23=R23*H(2,:)'/norm(H(2,:));
 F12=R12*H(3,:)'/norm(H(3,:));
 F34=R34*H(4,:)'/norm(H(4,:)); 
 
 F=F14+F23+F12+F34;
 
 %% Ps matrix --- control function base matrix
omega=2*pi/T;
P_vec1=[1 sin(omega*t) cos(omega*t) sin(omega*2*t) cos(omega*2*t)];
P_s=blkdiag(P_vec1, P_vec1);

 %% control system
 B=[zeros(2, 3) eye(2)]';
g=[zeros(5,2); (1/R)*P\B];
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
f=[x_two; P\(-D+F)];
A_lin=jacobian(f+ g*P_s*lam, x);
matlabFunction(f, 'file', 'sfun_f', 'vars', {x, epstau});
matlabFunction(A_lin, 'file', 'sfun_A', 'vars', {x, epstau});