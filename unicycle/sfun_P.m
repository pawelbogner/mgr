function P_s = sfun_P(t)
%SFUN_P
%    P_S = SFUN_P(T)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    08-Apr-2014 11:46:21

t2 = pi.*t.*(1.0./2.0);
t3 = sin(t2);
t4 = cos(t2);
P_s = reshape([1.0,0.0,t3,0.0,t4,0.0,0.0,1.0,0.0,t3,0.0,t4],[2, 6]);