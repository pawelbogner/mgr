function A = mp_inverse(P)
    PP=P*P';
     deter=det(PP)
     if(deter<10e-6)
         PP=PP+0.0001*eye(size(PP));
         deter2=det(PP)
    A= P'/(PP);
end