close all

%% Dinamica, Long Ambiguity monosettoriale

periods=    	100000;

% T=periods*10; % use T to find numerically steady states

dt=1/1000;


phi=1.22;

psi=1;
rho=.03;
delta=.1;
al=1/3;
lambda=5;
eta=200;
A=1;

omega=al*(1+phi)/(al+phi);
mu=omega/al;
Theta=(eta/(eta-1))^al*lambda^(al/eta);

gamma=al/(al-eta*(al-1));
Omega=(al/rho)^(gamma*(eta-1))*( A^(eta/al)*(eta/(eta-1))^(1/al)*lambda     )^(gamma);

x1=     1;
% x=      [x1;zeros(T-1,1)];
% x=x(:);
% b=x*0;
% k=x+b;
% y=x;
% r=x;
% l=x;
% w=x;
% zbar=x;
% c=x;

%% zbar and l are time invariant (also conditional of TFP or broad productivity)

zbar_TFP_REP=A*0+lambda^(inv(eta));
l_TFP_REP=A*0+inv(psi)^inv(1+phi);

zbar_PPP=zbar_TFP_REP;
l_PPP=l_TFP_REP;

%% Dynamics standard TFP, use for SS

% for n=1:T
%     
%     Theta= (eta*inv(eta-1))^al*lambda^(al/eta);
% 
%     y(n)=Theta*A*x(n)^al*l(n)^(1-al);
%     
%     r(n)=al*(eta-1)/eta*y(n)/x(n);
% 
%     w(n)=(1-al)*(n)/l(n);
% 
%     c(n)=w(n)*l(n);
%     
%     xdot=al/eta*y(n)+(r(n)-delta)*x(n);
%     
%     x(n+1)=x(n)+xdot*dt;
%     
% end


x_TFP_REP=((eta*inv(eta-1))^al*lambda^(al/eta)*A*al/delta)^inv(1-al)*l_TFP_REP(1);

%% Take SS

x_TFP_REP=      [x_TFP_REP(end);zeros(periods-1,1)];
x_TFP_REP=x_TFP_REP(:);
b_TFP_REP=x_TFP_REP*0;
k_TFP_REP=x_TFP_REP+b_TFP_REP;
y_TFP_REP=x_TFP_REP;
r_TFP_REP=x_TFP_REP;
l_TFP_REP=x_TFP_REP;
w_TFP_REP=x_TFP_REP;
zbar_TFP_REP=x_TFP_REP;
c_TFP_REP=x_TFP_REP;

x_PPP=x_TFP_REP(:);
b_PPP=x_TFP_REP*0;
k_PPP=x_PPP+b_PPP;
y_PPP=x_TFP_REP;
r_PPP=x_TFP_REP;
l_PPP=x_TFP_REP;
w_PPP=x_TFP_REP;
zbar_PPP=x_TFP_REP;
c_PPP=x_TFP_REP;

zbar_TFP_REP=zbar_TFP_REP*0+lambda^(inv(eta));
l_TFP_REP=l_TFP_REP*0+inv(psi)^inv(1+phi);

zbar_PPP=zbar_TFP_REP;
l_PPP=l_TFP_REP;

%% IR

P=A+x_TFP_REP*0;   %P=productivity

P(floor(periods/10))=1.1*P(floor(periods/10));  % 10percent IR

for n=(floor(periods/10)):periods

    if n==1
        Pdot=0;
    else
    Pdot=.10*(P(1)-P(n));
    end
    P(n+1)=P(n)+Pdot*dt;

end

AIR=P;

betagamma=1;

centering = -log(2)-betagamma*P(1);  % serve per far venire, in media al=1/3

alIR= exp(betagamma*P+centering)./(1+exp(betagamma*P+centering));  % alpha depends on P, productivity, through a mapping bounded between 0 and 1

%% IRF standard TFP

for n=1:periods

    A=AIR(n);

    Theta= (eta*inv(eta-1))^al*lambda^(al/eta);

    y_TFP_REP(n)=Theta*A*x_TFP_REP(n)^al*l_TFP_REP(n)^(1-al);
    
    r_TFP_REP(n)=al*(eta-1)/eta*y_TFP_REP(n)/x_TFP_REP(n);

    w_TFP_REP(n)=(1-al)*y_TFP_REP(n)/l_TFP_REP(n);

    c_TFP_REP(n)=w_TFP_REP(n)*l_TFP_REP(n);
    
    xdot=al/eta*y_TFP_REP(n)+(r_TFP_REP(n)-delta)*x_TFP_REP(n);
    
    x_TFP_REP(n+1)=x_TFP_REP(n)+xdot*dt;

end


%% IRF Productivity, P

for n=1:periods

    A=AIR(n);

    al=alIR(n);

    Theta= (eta*inv(eta-1))^al*lambda^(al/eta);

    y_PPP(n)=Theta*A*x_PPP(n)^al*l_PPP(n)^(1-al);
    
    r_PPP(n)=al*(eta-1)/eta*y_PPP(n)/x_PPP(n);

    w_PPP(n)=(1-al)*y_PPP(n)/l_PPP(n);

    c_PPP(n)=w_PPP(n)*l_PPP(n);
    
    xdot_PPP=al/eta*y_PPP(n)+(r_PPP(n)-delta)*x_PPP(n);
    
    x_PPP(n+1)=x_PPP(n)+xdot_PPP*dt;

end

%% IR definizione shock

PIR=(P-P(1))/P(1)*100;
AIR=(AIR-AIR(1))/AIR(1)*100; %classico
alIR=(alIR-alIR(1))/alIR(1)*100;


alIR_PPP=alIR;
alIR=alIR*0;

%% IR definizione TFP

xIR_TFP_REP=(x_TFP_REP-x_TFP_REP(1))/x_TFP_REP(1)*100;
yIR_TFP_REP=(y_TFP_REP-y_TFP_REP(1))/y_TFP_REP(1)*100;
cIR_TFP_REP=(c_TFP_REP-c_TFP_REP(1))/c_TFP_REP(1)*100;
zbarIR_TFP_REP=(zbar_TFP_REP-zbar_TFP_REP(1))/zbar_TFP_REP(1)*100;
lIR_TFP_REP=(l_TFP_REP-l_TFP_REP(1))/l_TFP_REP(1)*100;
wIR_TFP_REP=(w_TFP_REP-w_TFP_REP(1))/w_TFP_REP(1)*100;
rIR_TFP_REP=(r_TFP_REP-r_TFP_REP(1))/r_TFP_REP(1)*100;

xIR_TFP_REP=xIR_TFP_REP(1:end-1);

%% IR definizione P

xIR_PPP=(x_PPP-x_PPP(1))/x_PPP(1)*100;
yIR_PPP=(y_PPP-y_PPP(1))/y_PPP(1)*100;
cIR_PPP=(c_PPP-c_PPP(1))/c_PPP(1)*100;
zbarIR_PPP=(zbar_PPP-zbar_PPP(1))/zbar_PPP(1)*100;
lIR_PPP=(l_PPP-l_PPP(1))/l_PPP(1)*100;
wIR_PPP=(w_PPP-w_PPP(1))/w_PPP(1)*100;
rIR_PPP=(r_PPP-r_PPP(1))/r_PPP(1)*100;

xIR_PPP=xIR_PPP(1:end-1);



Points=60;

discretizzazione=linspace(1,periods,Points);
discretizzazione=round(discretizzazione);

k_TFP_REP=k_TFP_REP(discretizzazione);
x_TFP_REP=x_TFP_REP(discretizzazione);
b_TFP_REP=b_TFP_REP(discretizzazione);
y_TFP_REP=y_TFP_REP(discretizzazione);
r_TFP_REP=r_TFP_REP(discretizzazione);
l_TFP_REP=l_TFP_REP(discretizzazione);
w_TFP_REP=w_TFP_REP(discretizzazione);
zbar_TFP_REP=zbar_TFP_REP(discretizzazione);
c_TFP_REP=c_TFP_REP(discretizzazione);

xIR_PPP=xIR_PPP(discretizzazione);
yIR_PPP=yIR_PPP(discretizzazione);
cIR_PPP=cIR_PPP(discretizzazione);
wIR_PPP=wIR_PPP(discretizzazione);
rIR_PPP=rIR_PPP(discretizzazione);

xIR_TFP_REP=xIR_TFP_REP(discretizzazione);
yIR_TFP_REP=yIR_TFP_REP(discretizzazione);
cIR_TFP_REP=cIR_TFP_REP(discretizzazione);
wIR_TFP_REP=wIR_TFP_REP(discretizzazione);
rIR_TFP_REP=rIR_TFP_REP(discretizzazione);
