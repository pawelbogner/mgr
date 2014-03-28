function [s14, s23, s12, s34] = calc_slips(x)
N=size(x, 1);
s14=zeros(1,N);
s23=zeros(1,N);
s12=zeros(1,N);
s34=zeros(1,N);
a=0.730;
b=0.350;
for i=1:N
 H = [-sin(x(i,3)/a) cos(x(i,3)/a)    0  0   0 ;
      -sin(x(i,3)/a) cos(x(i,3)/a)    1  0   0 ;
      -sin(x(i,3)/a) cos(x(i,3)/a) -b/a -1   0 ;
      -sin(x(i,3)/a) cos(x(i,3)/a) -b/a  0  -1 ];
 
  %slips
  s14(i)=H(1,:)*x(i,6:10)';
  s23(i)=H(2,:)*x(i,6:10)';
  s12(i)=H(3,:)*x(i,6:10)';
  s34(i)=H(4,:)*x(i,6:10)';
end