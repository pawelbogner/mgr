function A_lin2 = sfun_A2(in1,in2)
%SFUN_A2
%    A_LIN2 = SFUN_A2(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    09-Mar-2014 17:44:16

q3 = in1(3,:);
qd1 = in1(5,:);
qd2 = in1(6,:);
qd4 = in1(8,:);
t4 = cos(q3);
t2 = abs(t4);
t6 = sin(q3);
t3 = abs(t6);
t5 = t2.^2;
t7 = t3.^2;
t8 = t5+t7;
t9 = conj(q3);
t10 = sqrt(t8);
t11 = cos(t9);
t12 = t5+t7+9.0./1.0e2;
t13 = sin(t9);
t14 = sqrt(t12);
t15 = 1.0./sqrt(t8);
t16 = sign(t4);
t17 = sign(t6);
t19 = t2.*t6.*t16.*2.0;
t20 = t3.*t4.*t17.*2.0;
t18 = t19-t20;
t21 = 1.0./sqrt(t12);
t22 = qd1.*t4.*t10.*t11.*1.0e1;
t23 = qd2.*t6.*t10.*t11.*1.0e1;
t24 = qd2.*t4.*t13.*t14.*1.0e4;
t25 = qd1.*t6.*t13.*t14.*1.0e4;
t27 = qd4.*t10.*t11.*3.0;
t26 = t22+t23+t24+t25-t27;
t28 = 1.0./t12.^(3.0./2.0);
t29 = qd1.*t4.*t10.*t13.*1.0e1;
t30 = qd2.*t6.*t10.*t13.*1.0e1;
t31 = 1.0./t8.^(3.0./2.0);
t32 = qd2.*t4.*t11.*t14.*1.0e4;
t33 = qd1.*t6.*t11.*t14.*1.0e4;
t34 = t29+t30+t32+t33-qd4.*t10.*t13.*3.0;
A_lin2 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t15.*t21.*(t29+t30-qd4.*t10.*t13.*3.0-qd2.*t4.*t10.*t11.*1.0e1+qd1.*t6.*t10.*t11.*1.0e1-qd2.*t4.*t11.*t14.*1.0e4-qd1.*t4.*t13.*t14.*1.0e4-qd1.*t6.*t11.*t14.*1.0e4+qd2.*t6.*t13.*t14.*1.0e4-qd4.*t11.*t15.*t18.*(3.0./2.0)+qd1.*t4.*t11.*t15.*t18.*5.0+qd2.*t6.*t11.*t15.*t18.*5.0+qd2.*t4.*t13.*t18.*t21.*5.0e3+qd1.*t6.*t13.*t18.*t21.*5.0e3).*2.943e-2-t15.*t18.*t26.*t28.*1.4715e-2-t18.*t21.*t26.*t31.*1.4715e-2,t15.*t21.*(-t22-t23+t24+t25+t27-qd2.*t4.*t10.*t13.*1.0e1-qd1.*t4.*t11.*t14.*1.0e4+qd1.*t6.*t10.*t13.*1.0e1+qd2.*t6.*t11.*t14.*1.0e4-qd4.*t13.*t15.*t18.*(3.0./2.0)+qd1.*t4.*t13.*t15.*t18.*5.0+qd2.*t6.*t13.*t15.*t18.*5.0+qd2.*t4.*t11.*t18.*t21.*5.0e3+qd1.*t6.*t11.*t18.*t21.*5.0e3).*2.943e-2-t15.*t18.*t28.*t34.*1.4715e-2-t18.*t21.*t31.*t34.*1.4715e-2,0.0,t21.*(qd2.*t4.*1.0e1-qd1.*t6.*1.0e1).*1.962e-1+t18.*t28.*(qd4.*-3.0+qd1.*t4.*1.0e1+qd2.*t6.*1.0e1).*9.81e-2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,t15.*t21.*(t4.*t10.*t11.*1.0e1+t6.*t13.*t14.*1.0e4).*(-2.943e-2),t15.*t21.*(t4.*t10.*t13.*1.0e1+t6.*t11.*t14.*1.0e4).*(-2.943e-2),0.0,t4.*t21.*(9.81e2./5.0e2),0.0,1.0,0.0,0.0,t15.*t21.*(t6.*t10.*t11.*1.0e1+t4.*t13.*t14.*1.0e4).*(-2.943e-2),t15.*t21.*(t4.*t11.*t14.*1.0e4+t6.*t10.*t13.*1.0e1).*(-2.943e-2),0.0,t6.*t21.*(9.81e2./5.0e2),0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,t11.*t21.*8.829e-2,t13.*t21.*8.829e-2,0.0,t21.*(-5.886e-1)],[8, 8]);
