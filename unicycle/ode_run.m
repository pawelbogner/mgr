function [tout, xout] = ode_run(lambdas)

opt = odeset('Events', @event);

tstart=0;
tfinal=4;
x0=[0; 0; pi/2; zeros(5,1); zeros(8*2*3,1)];

tout=tstart;
xout=x0';
teout=[];
xeout=[];
ieout=[];

skids=zeros(2,1);
counter=0;
ffcn=@sfun_f0;
Afcn=@sfun_A0;
while tout(end) < tfinal
    % ode_fun = [f+g Ps lam; A ksi+ BP]

    odefcn=@(t,x) [ffcn(x(1:8))+sfun_g(x(1:8))*sfun_P(t)*lambdas; reshape(Afcn(x(1:8), lambdas)*reshape(x(9:end), 8, [])+sfun_B(x(1:8))*sfun_P(t), [], 1)];
    
    %f= @(t,x) ode_function(t, x, skids, lambdas);
    [t,x,te,xe,ie] = ode45(odefcn, [tstart, tfinal], x0, opt);
    %ie
    tlen=length(t);
    tout=[tout; t(2:end)];
    xout=[xout; x(2:end,:)];
    teout=[teout; te];
    xeout=[xeout; xe];
    ieout=[ieout; ie];
    
    x0=x(end,:);
    tstart=t(end);
    
    %tutaj sprawdzamy gdzie zmienialy sie poslizgi
    for k=1:length(ie)
        if ie(k)<=2 %skid start
            skids(ie(k))=1;
        else %skid end
            skids(ie(k)-2)=0;
        end
        
    end
    %wyliczamy aktualny index poslizgow
    skid_index=0;
    for k=0:1
        skid_index=skid_index+bitshift(1, k)*skids(1+k);
    end
        
    %tutaj robimy przelaczanie
    ffcn=str2func(['sfun_f' int2str(skid_index)]);
    Afcn=str2func(['sfun_A' int2str(skid_index)]);
    counter=counter+1;
    
    %disp('bylo przelaczenie');
    %[counter tstart]
    %plot(xout(:,1), xout(:,2))
    
end


