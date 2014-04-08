function s = calc_slips(x)
N=size(x, 1);
s1=zeros(1,N);
s2=zeros(1,N);
R=0.3;
for i=1:N
  H=[sin(x(i,3)) -cos(x(i,3)) 0  0;
     cos(x(i,3))  sin(x(i,3)) 0 -R];
 
  %slips
  s1(i)=H(1,:)*x(i,5:8)';
  s2(i)=H(2,:)*x(i,5:8)';
end

s=[s1', s2'];