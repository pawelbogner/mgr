function k=kin(q)
a2=0.1276;
a3=0.041;
a4=0.2615;
a5=0.0195;

a23=a2+a3;
a4a5c5=a4+a5*cos(q(5));
q24=q(2)+q(4);
c24=cos(q24);
a4a5c5c24=a4a5c5*c24;
c2=cos(q(2));
a23c2=a23*c2;
a4a5c5c24a23c2=a4a5c5c24+a23c2;
s5=sin(q(5));
s1=sin(q(1));
c1=sin(q(1));
a5s5=a5*s5;

x=-a5s5*s1+ c1*(a4a5c5c24a23c2);
y= a5s5*c1+ s1*(a4a5c5c24a23c2);
z= a23*sin(q(2))+ a4a5c5*sin(q24);
phi=q(5)+pi/2;
theta=q24;
psi=q(1)-pi/2;

k=[x; y; z; phi; theta; psi];