function A = mp_inverse(P)
    PP=P*P.';
     deter=det(PP)
%       if(deter<10e-6)
%           PP=PP+0.005*eye(size(PP));
%           deter2=det(PP)
%       end
    A= P.'/(PP);
end