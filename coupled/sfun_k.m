function k = sfun_k(in1)
%SFUN_K
%    K = SFUN_K(IN1)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    26-Mar-2014 20:33:11

x_one1 = in1(1,:);
x_one2 = in1(2,:);
x_one3 = in1(3,:);
k = [x_one1;x_one2;x_one3.*(1.0e2./7.3e1)];
