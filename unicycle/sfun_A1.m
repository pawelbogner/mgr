function A_lin1 = sfun_A1(in1,in2)
%SFUN_A1
%    A_LIN1 = SFUN_A1(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    08-Apr-2014 11:46:24

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
t22 = qd1.*t4.*t10.*t11.*3.0e1;
t23 = qd2.*t6.*t10.*t11.*3.0e1;
t24 = qd1.*t6.*t13.*t14;
t26 = qd4.*t10.*t11.*9.0;
t27 = qd2.*t4.*t13.*t14;
t25 = t22+t23+t24-t26-t27;
t28 = 1.0./t12.^(3.0./2.0);
t29 = qd4.*t10.*t13.*9.0;
t30 = qd1.*t6.*t11.*t14;
t31 = 1.0./t8.^(3.0./2.0);
t32 = qd2.*t4.*t11.*t14;
t33 = qd1.*t4.*t10.*t13.*3.0e1;
t34 = qd2.*t6.*t10.*t13.*3.0e1;
t35 = -t29-t30+t32+t33+t34;
A_lin1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t15.*t21.*(t29+t30+qd2.*t4.*t10.*t11.*3.0e1-qd1.*t4.*t10.*t13.*3.0e1-qd1.*t6.*t10.*t11.*3.0e1-qd2.*t4.*t11.*t14-qd2.*t6.*t10.*t13.*3.0e1+qd1.*t4.*t13.*t14+qd2.*t6.*t13.*t14+qd4.*t11.*t15.*t18.*(9.0./2.0)-qd1.*t4.*t11.*t15.*t18.*1.5e1-qd2.*t6.*t11.*t15.*t18.*1.5e1+qd2.*t4.*t13.*t18.*t21.*(1.0./2.0)-qd1.*t6.*t13.*t18.*t21.*(1.0./2.0)).*(-9.81e2./1.0e3)-t15.*t18.*t25.*t28.*4.905e-1-t18.*t21.*t25.*t31.*4.905e-1,t15.*t21.*(-t22-t23-t24+t26+t27-qd2.*t4.*t10.*t13.*3.0e1+qd1.*t4.*t11.*t14+qd1.*t6.*t10.*t13.*3.0e1+qd2.*t6.*t11.*t14-qd4.*t13.*t15.*t18.*(9.0./2.0)+qd1.*t4.*t13.*t15.*t18.*1.5e1+qd2.*t6.*t13.*t15.*t18.*1.5e1+qd2.*t4.*t11.*t18.*t21.*(1.0./2.0)-qd1.*t6.*t11.*t18.*t21.*(1.0./2.0)).*(9.81e2./1.0e3)-t15.*t18.*t28.*t35.*4.905e-1-t18.*t21.*t31.*t35.*4.905e-1,0.0,t21.*(qd2.*t4.*1.0e1-qd1.*t6.*1.0e1).*(9.81e2./5.0e1)+t18.*t28.*(qd4.*-3.0+qd1.*t4.*1.0e1+qd2.*t6.*1.0e1).*(9.81e2./1.0e2),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,t15.*t21.*(t4.*t10.*t11.*3.0e1+t6.*t13.*t14).*(-9.81e2./1.0e3),t15.*t21.*(t4.*t10.*t13.*3.0e1-t6.*t11.*t14).*(-9.81e2./1.0e3),0.0,t4.*t21.*(9.81e2./5.0),0.0,1.0,0.0,0.0,t15.*t21.*(t6.*t10.*t11.*3.0e1-t4.*t13.*t14).*(-9.81e2./1.0e3),t15.*t21.*(t4.*t11.*t14+t6.*t10.*t13.*3.0e1).*(-9.81e2./1.0e3),0.0,t6.*t21.*(9.81e2./5.0),0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,t11.*t21.*8.829,t13.*t21.*8.829,0.0,t21.*(-5.886e1)],[8, 8]);
