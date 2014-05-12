%% variables
n=5; %dim x
x=sym('x', [2*n 1]);
m=2;
%u=sym('u', [m 1]);
s=m*(2*4+1);
%s=m*4;

lam=sym('lam', [s 1]);
syms t;

 r=4;
 q=sym('q', [r 1]);

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
ap2=0.008;
gr=9.81;
T_h=20;

%% Q elements
Q11=mp+4*mw;
Q22=Q11;

Q44=2*Iw33;
Q55=Q44;

Q33= Ip33 + mp*(ap1.^2+ap2.^2) + 4*(Iw11 + mw*b.^2) + 2*mw*a.^2;

Q13 = -mp*(ap1*sin(x(3)./a) + ap2*cos(x(3)./a)) - 2*mw*a*sin(x(3)./a);
Q23 =  mp*(ap1*cos(x(3)./a) - ap2*sin(x(3)./a)) + 2*mw*a*cos(x(3)./a);

%% P matrix
P = [Q11    0      Q13./a      0           0;
     0      Q22    Q23./a      0           0;
     Q13./a Q23./a Q33./(a.^2) 0           0;
     0      0      0           Q44./(R.^2) 0;
     0      0      0           0           Q55./(R.^2)];
 
 %% D matrix
 D = ((x(8).^2)./(a.^2)).*([-Q23; Q13; 0; 0; 0]);
 
 %% H matrix
  
   H1=[-sin(x(3)./a) cos(x(3)./a)    0  0   0];
   H2=[-sin(x(3)./a) cos(x(3)./a)    1  0   0];
   H3=[ cos(x(3)./a) sin(x(3)./a) -b/a -1   0];
   H4=[ cos(x(3)./a) sin(x(3)./a)  b/a  0  -1];
   
   H1T=[-sin(x(3)./a); cos(x(3)./a);    0;  0;   0];
   H2T=[-sin(x(3)./a); cos(x(3)./a);    1;  0;   0];
   H3T=[ cos(x(3)./a); sin(x(3)./a); -b/a; -1;   0];
   H4T=[ cos(x(3)./a); sin(x(3)./a);  b/a;  0;  -1];
 
  %slips
  s14=H1*x(6:10);
  s23=H2*x(6:10);
  s12=H3*x(6:10);
  s34=H4*x(6:10);
  
 %% F vector
% N forces
 N1=((mp/4)+mw)*gr;
 N2=((mp/4)+mw)*gr;
 N3=((mp/4)+mw)*gr;
 N4=((mp/4)+mw)*gr;
 
% R coefficients
 R14=-(1*N1 + 1*N4)*s14;
 R23=-(1*N2 + 1*N3)*s23;
 R12=-(3*N1 + 3*N2)*s12;
 R34=-(3*N3 + 3*N4)*s34;
 
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
omega=2*3.141592653589793/T_h;
P_s=[1 sin(omega*t) cos(omega*t) sin(2*omega*t) cos(2*omega*t) sin(3*omega*t) cos(3*omega*t) sin(4*omega*t) cos(4*omega*t) 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 1 sin(omega*t) cos(omega*t) sin(2*omega*t) cos(2*omega*t) sin(3*omega*t) cos(3*omega*t) sin(4*omega*t) cos(4*omega*t)];
 
  %% Legendre
%  l0=1;
%  l1=2/T_h^2*(2*t - T_h);
%  l2=4/T_h^4*(6*t^2 - 6*t*T_h + T_h^2);
%  l3=8/T_h^6*(2*t - T_h)*(10*t^2 - 10*t*T_h + T_h^2);
%  
%  P_s=[l0 l1 l2 l3 0 0 0 0;
%       0 0 0 0 l0 l1 l2 l3];

 %% control system
 B=[0 0;
    0 0;
    0 0;
    1 0;
    0 1];
g=[zeros(5,2); inv(P)*B./R];
B_lin=g;

%% output function

a2=0.1276;
a3=0.041;
a4=0.2615;
a5=0.0195;

c1=cos(q(1));
s1=sin(q(1));
c2=cos(q(2));
s2=sin(q(2));
c24=cos(q(2)+q(3));
s24=sin(q(2)+q(3));
c5=cos(q(4));
s5=sin(q(4));

x_m=c1*((a2+a3)*c2 + c24*(a4+a5*c5)) - a5*s1*s5;
y_m=s1*((a2+a3)*c2 + c24*(a4+a5*c5)) + a5*c1*s5;
z_m=    (a2+a3)*s2 + s24*(a4+a5*c5);

phi=q(4)+pi/2;
theta=q(2)+q(3);
psi=q(1)-pi/2;

k=[x(1); x(2); x(3)./a; x(6); x(7); x(8)./a; z_m; phi; theta; psi];

%%
C_lin=jacobian(k, x);
D_lin=jacobian(k, q);
matlabFunction(k, 'file', 'sfun_k', 'vars', {x, q});
matlabFunction(B_lin, 'file', 'sfun_B', 'vars', {x});
matlabFunction(C_lin, 'file', 'sfun_C', 'vars', {x});
matlabFunction(D_lin, 'file', 'sfun_D', 'vars', {q});
matlabFunction(g, 'file', 'sfun_g', 'vars', {x});
matlabFunction(P_s, 'file', 'sfun_P', 'vars', {t});


 %% generation of A_lin and f
f=[x(6:10); inv(P)*(-D+F)];
A_lin=jacobian(f+ g*P_s*lam, x);
matlabFunction(f, 'file', 'sfun_f', 'vars', {x});
matlabFunction(A_lin, 'file', 'sfun_A', 'vars', {x});