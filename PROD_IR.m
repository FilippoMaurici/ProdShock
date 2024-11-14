clear all
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
eta=2;
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

zbar=A*0+lambda^(inv(eta));
l=A*0+inv(psi)^inv(1+phi);

zbar_prod=zbar;
l_prod=l;

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


x=((eta*inv(eta-1))^al*lambda^(al/eta)*A*al/delta)^inv(1-al)*l(1);

%% Take SS

x=      [x(end);zeros(periods-1,1)];
x=x(:);
b=x*0;
k=x+b;
y=x;
r=x;
l=x;
w=x;
zbar=x;
c=x;

x_prod=x(:);
b_prod=x*0;
k_prod=x_prod+b_prod;
y_prod=x;
r_prod=x;
l_prod=x;
w_prod=x;
zbar_prod=x;
c_prod=x;

zbar=zbar*0+lambda^(inv(eta));
l=l*0+inv(psi)^inv(1+phi);

zbar_prod=zbar;
l_prod=l;

%% IR

P=A+x*0;   %P=productivity

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

    y(n)=Theta*A*x(n)^al*l(n)^(1-al);
    
    r(n)=al*(eta-1)/eta*y(n)/x(n);

    w(n)=(1-al)*y(n)/l(n);

    c(n)=w(n)*l(n);
    
    xdot=al/eta*y(n)+(r(n)-delta)*x(n);
    
    x(n+1)=x(n)+xdot*dt;

end


%% IRF Productivity, P

for n=1:periods

    A=AIR(n);

    al=alIR(n);

    Theta= (eta*inv(eta-1))^al*lambda^(al/eta);

    y_prod(n)=Theta*A*x_prod(n)^al*l_prod(n)^(1-al);
    
    r_prod(n)=al*(eta-1)/eta*y_prod(n)/x_prod(n);

    w_prod(n)=(1-al)*y_prod(n)/l_prod(n);

    c_prod(n)=w_prod(n)*l_prod(n);
    
    xdot_prod=al/eta*y_prod(n)+(r_prod(n)-delta)*x_prod(n);
    
    x_prod(n+1)=x_prod(n)+xdot_prod*dt;

end

%% IR definizione shock

PIR=(P-P(1))/P(1)*100;
AIR=(AIR-AIR(1))/AIR(1)*100; %classico
alIR=(alIR-alIR(1))/alIR(1)*100;


alIR_prod=alIR;
alIR=alIR*0;

%% IR definizione TFP

xIR=(x-x(1))/x(1)*100;
yIR=(y-y(1))/y(1)*100;
cIR=(c-c(1))/c(1)*100;
zbarIR=(zbar-zbar(1))/zbar(1)*100;
lIR=(l-l(1))/l(1)*100;
wIR=(w-w(1))/w(1)*100;
rIR=(r-r(1))/r(1)*100;

xIR=xIR(1:end-1);

%% IR definizione P

xIR_prod=(x_prod-x_prod(1))/x_prod(1)*100;
yIR_prod=(y_prod-y_prod(1))/y_prod(1)*100;
cIR_prod=(c_prod-c_prod(1))/c_prod(1)*100;
zbarIR_prod=(zbar_prod-zbar_prod(1))/zbar_prod(1)*100;
lIR_prod=(l_prod-l_prod(1))/l_prod(1)*100;
wIR_prod=(w_prod-w_prod(1))/w_prod(1)*100;
rIR_prod=(r_prod-r_prod(1))/r_prod(1)*100;

xIR_prod=xIR_prod(1:end-1);



%% Plots

Points=60;

discretizzazione=linspace(1,periods,Points);
discretizzazione=round(discretizzazione);

k=k(discretizzazione);
x=x(discretizzazione);
b=b(discretizzazione);
y=y(discretizzazione);
r=r(discretizzazione);
l=l(discretizzazione);
w=w(discretizzazione);
zbar=zbar(discretizzazione);
c=c(discretizzazione);

xIR_prod=xIR_prod(discretizzazione);
yIR_prod=yIR_prod(discretizzazione);
cIR_prod=cIR_prod(discretizzazione);
wIR_prod=wIR_prod(discretizzazione);
rIR_prod=rIR_prod(discretizzazione);

xIR=xIR(discretizzazione);
yIR=yIR(discretizzazione);
cIR=cIR(discretizzazione);
wIR=wIR(discretizzazione);
rIR=rIR(discretizzazione);

PIR=PIR(discretizzazione);
AIR=AIR(discretizzazione);
alIR_prod=alIR_prod(discretizzazione);
alIR=alIR(discretizzazione);

ylim padded

close all

set(gca,'Color','none')

g1 = figure('Name','PIR');
 
 plot(discretizzazione,PIR,'r--','linewidth' , 2)
 hold on
 plot(discretizzazione,PIR,'b-o','linewidth' , 1)

 
% plot(discretizzazione,CUBueraIR,'k-.','linewidth' , 1)
h = yline(0, 'k', 'LineWidth', 1, 'HandleVisibility','off');
legend('$\mathbf{P}$, TFP shock','$\mathbf{P}$, productivity shock','Interpreter','latex','Location','best')
%xlabel('$t$','Interpreter','latex')
%%title('Producing entrepreneur consumption, IRF','Interpreter','latex')
xlim([0 max(discretizzazione)]) %([0 size(PIR,1)])
%%ylim([1.1*min(xIR)-.1 1.1*max(xIR)+.1])
set(gca,'Color','none')
% set(gca,'XTick',[])
ylim padded

print(g1,'PIR','-depsc',  '-painters','-r600')
% export_fig('avgCUIRF', '-dpngc', '-transparent', '-r300');

close all

set(gca,'Color','none')
g1 = figure('Name','AIR');



plot(discretizzazione,AIR,'r--','linewidth' , 2)
hold on
plot(discretizzazione,AIR,'b-o','linewidth' , 1)


% plot(discretizzazione,CUBueraIR,'k-.','linewidth' , 1)
h = yline(0, 'k', 'LineWidth', 1, 'HandleVisibility','off');
legend('$A$, TFP shock','$A$, productivity shock','Interpreter','latex','Location','best')
%xlabel('$t$','Interpreter','latex')
%%title('Producing entrepreneur consumption, IRF','Interpreter','latex')
xlim([0 max(discretizzazione)]) %([0 max(discretizzazione)]) %([0 size(xIR,1)])
%%ylim([1.1*min(xIR)-.1 1.1*max(xIR)+.1])
set(gca,'Color','none')
% set(gca,'XTick',[])
ylim padded

print(g1,'AIR','-depsc',  '-painters','-r600')
% export_fig('avgCUIRF', '-dpngc', '-transparent', '-r300');

close all

set(gca,'Color','none')
g1 = figure('Name','alphaIR');


plot(discretizzazione,alIR,'r--','linewidth' , 2)
hold on
plot(discretizzazione,alIR_prod,'b-o','linewidth' , 1)


% plot(discretizzazione,CUBueraIR,'k-.','linewidth' , 1)
hold on
h = yline(0, 'k', 'LineWidth', 1, 'HandleVisibility','off');
legend('$\alpha$, TFP shock','$\alpha$, productivity shock','Interpreter','latex','Location','best')
%xlabel('$t$','Interpreter','latex')
%%title('Producing entrepreneur consumption, IRF','Interpreter','latex')
xlim([0 max(discretizzazione)]) %([0 size(xIR,1)])
%%ylim([1.1*min(xIR)-.1 1.1*max(xIR)+.1])
set(gca,'Color','none')
% set(gca,'XTick',[])
ylim padded

print(g1,'alphaIR','-depsc',  '-painters','-r600')
% export_fig('avgCUIRF', '-dpngc', '-transparent', '-r300');


close all

set(gca,'Color','none')
g1 = figure('Name','yIR');



plot(discretizzazione,yIR,'r--','linewidth' , 2)
hold on
plot(discretizzazione,yIR_prod,'b-o','linewidth' , 1)


% plot(discretizzazione,CUBueraIR,'k-.','linewidth' , 1)
h = yline(0, 'k', 'LineWidth', 1, 'HandleVisibility','off');
legend('$y$, TFP shock','$y$, productivity shock','Interpreter','latex','Location','best')

xlim([0 max(discretizzazione)]) %([0 size(yIR,1)])
%%ylim([1.1*min(yIR)-.1 1.1*max(yIR)+.1])
set(gca,'Color','none')
% set(gca,'XTick',[])
ylim padded

print(g1,'yIR','-depsc',  '-painters','-r600')
% export_fig('avgCUIRF', '-dpngc', '-transparent', '-r300');

close all

set(gca,'Color','none')
g1 = figure('Name','wIR');


plot(discretizzazione,wIR,'r--','linewidth' , 2)
hold on
plot(discretizzazione,wIR_prod,'b-o','linewidth' , 1)

 hold on

% plot(discretizzazione,CUBueraIR,'k-.','linewidth' , 1)
h = yline(0, 'k', 'LineWidth', 1, 'HandleVisibility','off');
legend('$w$, TFP shock','$w$, productivity shock','Interpreter','latex','Location','best')

%xlabel('$t$','Interpreter','latex')
%%title('Producing entrepreneur consumption, IRF','Interpreter','latex')
xlim([0 max(discretizzazione)]) %([0 size(wIR,1)])
%%ylim([1.1*min(wIR)-.1 1.1*max(wIR)+.1])
set(gca,'Color','none')
% set(gca,'XTick',[])
ylim padded

print(g1,'wIR','-depsc',  '-painters','-r600')
% export_fig('avgCUIRF', '-dpngc', '-transparent', '-r300');

close all

set(gca,'Color','none')
g1 = figure('Name','rIR');


plot(discretizzazione,rIR,'r--','linewidth' , 2)
hold on

plot(discretizzazione,rIR_prod,'b-o','linewidth' , 1)


% plot(discretizzazione,CUBueraIR,'k-.','linewidth' , 1)
h = yline(0, 'k', 'LineWidth', 1, 'HandleVisibility','off');
legend('$r$, TFP shock','$r$, productivity shock','Interpreter','latex','Location','best')

xlim([0 max(discretizzazione)]) %([0 size(rIR,1)])
%%ylim([1.1*min(rIR)-.1 1.1*max(rIR)+.1])
set(gca,'Color','none')
% set(gca,'XTick',[])
ylim padded

print(g1,'rIR','-depsc',  '-painters','-r600')
% export_fig('avgCUIRF', '-dpngc', '-transparent', '-r300');

close all

set(gca,'Color','none')
g1 = figure('Name','cIR');



plot(discretizzazione,cIR,'r--','linewidth' , 2)
hold on
plot(discretizzazione,cIR_prod,'b-o','linewidth' , 1)

 
% plot(discretizzazione,CUBueraIR,'k-.','linewidth' , 1)
h = yline(0, 'k', 'LineWidth', 1, 'HandleVisibility','off');
legend('$c$, TFP shock','$c$, productivity shock','Interpreter','latex','Location','best')

xlim([0 max(discretizzazione)]) %([0 size(cIR,1)])
%%ylim([1.1*min(cIR)-.1 1.1*max(cIR)+.1])
set(gca,'Color','none')
% set(gca,'XTick',[])
ylim padded

print(g1,'cIR','-depsc',  '-painters','-r600')
% export_fig('avgCUIRF', '-dpngc', '-transparent', '-r300');

close all

set(gca,'Color','none')
g1 = figure('Name','xIR');


plot(discretizzazione,xIR,'r--','linewidth' , 2) 
hold on
plot(discretizzazione,xIR_prod,'b-o','linewidth' , 1)


% plot(discretizzazione,CUBueraIR,'k-.','linewidth' , 1)
h = yline(0, 'k', 'LineWidth', 1, 'HandleVisibility','off');
legend('$x$, TFP shock','$x$, productivity shock','Interpreter','latex','Location','best')

xlim([0 max(discretizzazione)]) %([0 size(xIR,1)])
%%ylim([1.1*min(bIR)-.1 1.1*max(bIR)+.1])
set(gca,'Color','none')
% set(gca,'XTick',[])
ylim padded

print(g1,'xIR','-depsc',  '-painters','-r600')
