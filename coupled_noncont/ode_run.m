function [tout, xout] = ode_run(lambdas)

opt = odeset('Events', @event);

tstart=0;
tfinal=4;
%x0=[0; 0; 0.730*pi/2; zeros(7,1)]; %% dim Ksi= n*m*s
x0=[0; 0; 0.730*pi/2; zeros(7,1); zeros(10*2*7,1)];

tout=tstart;
xout=x0';
teout=[];
xeout=[];
ieout=[];

counter=0;
epstau=[30 30 30 30 30 30]';
while tout(end) < tfinal
    % ode_fun = [f+g Ps lam; A ksi+ BP]

    odefcn=@(t,x) [sfun_f(x(1:10), epstau)+sfun_g(x(1:10))*sfun_P(t)*lambdas; reshape(sfun_A(x(1:10), epstau)*reshape(x(11:end), 10, [])+sfun_B(x(1:10))*sfun_P(t), [], 1)];
    
    %f= @(t,x) ode_function(t, x, skids, lambdas);
    [t,x,te,xe,ie] = ode45(odefcn, [tstart, tfinal], x0, opt);
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
        if ie(k)<=4 %skid start
            epstau(ie(k))=30;
        else %skid end
            epstau(ie(k)-4)=0.1;
        end
    end
    %tutaj robimy przelaczanie
    
    counter=counter+1;
    
    disp('bylo przelaczenie');
    [counter tstart]
    %plot(xout(:,1), xout(:,2));
    
end
