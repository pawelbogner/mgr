function A = mp_inverse(P)
    PP=P*P';
     %deter=det(PP)
%     if(deter<10e-6)
%         PP=PP+1*eye(size(PP));
    A= P'/(PP);
end