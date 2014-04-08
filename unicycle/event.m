function [value,isterminal,direction] = event(~,q)
R=0.3;
H=[sin(q(3)) -cos(q(3)) 0  0;
   cos(q(3))  sin(q(3)) 0 -R];
 
%skids
s1=H(1,:)*q(5:8);
s2=H(2,:)*q(5:8);
s=[s1; s2];
% Locate the time when height passes through zero in a 
% decreasing direction and stop integration.
l=0.15;
value = [abs(s)-l; abs(s)-l];
isterminal = ones(2*length(s), 1);   % Stop the integration
direction = [ones(length(s), 1); -ones(length(s), 1)];   % Negative direction only