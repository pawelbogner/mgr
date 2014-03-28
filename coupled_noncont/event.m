function [value,isterminal,direction] = event(~,x)
a=0.730;
b=0.350;

 H = [-sin(x(3)/a) cos(x(3)/a)    0  0   0;
      -sin(x(3)/a) cos(x(3)/a)    1  0   0;
       cos(x(3)/a) sin(x(3)/a) -b/a -1   0;
       cos(x(3)/a) sin(x(3)/a)  b/a  0  -1];
 
  %skids
  s14=H(1,:)*x(6:10);
  s23=H(2,:)*x(6:10);
  s12=H(3,:)*x(6:10);
  s34=H(4,:)*x(6:10);
  s=[s14; s23; s12; s34];
% Locate the time when height passes through zero in a 
% decreasing direction and stop integration.
l=0.1;
value = [abs(s)-l; abs(s)-l];
isterminal = ones(2*length(s), 1);   % Stop the integration
direction = [ones(length(s), 1); -ones(length(s), 1)];   % Negative direction only