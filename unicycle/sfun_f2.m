function f2 = sfun_f2(in1)
%SFUN_F2
%    F2 = SFUN_F2(IN1)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    08-Apr-2014 11:46:20

q3 = in1(3,:);
qd1 = in1(5,:);
qd2 = in1(6,:);
qd3 = in1(7,:);
qd4 = in1(8,:);
t4 = cos(q3);
t2 = abs(t4);
t6 = sin(q3);
t3 = abs(t6);
t5 = t2.^2;
t7 = t3.^2;
t8 = t5+t7;
t9 = conj(q3);
t10 = cos(t9);
t11 = sqrt(t8);
t12 = t5+t7+9.0./1.0e2;
t13 = sin(t9);
t14 = sqrt(t12);
t15 = 1.0./sqrt(t8);
t16 = 1.0./sqrt(t12);
f2 = [qd1;qd2;qd3;qd4;t15.*t16.*(qd4.*t10.*t11.*-9.0+qd1.*t4.*t10.*t11.*3.0e1+qd2.*t6.*t10.*t11.*3.0e1-qd2.*t4.*t13.*t14.*1.0e2+qd1.*t6.*t13.*t14.*1.0e2).*(-9.81e-2);t15.*t16.*(qd4.*t11.*t13.*-9.0+qd1.*t4.*t11.*t13.*3.0e1+qd2.*t4.*t10.*t14.*1.0e2-qd1.*t6.*t10.*t14.*1.0e2+qd2.*t6.*t11.*t13.*3.0e1).*(-9.81e-2);0.0;t16.*(qd4.*-3.0+qd1.*t4.*1.0e1+qd2.*t6.*1.0e1).*(9.81e2./5.0e2)];
