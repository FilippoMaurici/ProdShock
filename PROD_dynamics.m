clear all
close all

%% Dinamica, Long Ambiguity monosettoriale

periods=    	100000;

dt=1/1000;

phi=1.22;
beta=1;
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

x1=     1;
x=      [x1;zeros(periods-1,1)];
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
% g=x;

%% Dynamics

zbar=zbar+lambda^(-eta);



for n=1:periods;


    
    w(n)=(1-al)/(1)*y(n)/l(n);
    r(n)=al*(eta-1)/((1)*eta)*y(n)/k(n);
    %     g(n)=w(n)*l(n)+r(n)*k(n);
    bdot=w(n)*l(n)+r(n)*b(n)-c(n);
    xdot=al/eta*y(n)+(r(n)-delta)*x(n);
    b(n+1)=b(n)+bdot*dt;
    x(n+1)=x(n)+xdot*dt;
    k(n+1)=x(n+1)+b(n+1);
end


for n=1:periods;

    Hx=(x_prod(n)-xSS_prod)*H(1);
    Hk=(k_prod(n)-kSS_prod)*H(2);

    c_prod(n)=-(Hx+Hk)/H(3)+cSS_prod;

    zbar_prod(n)=(lambda*x_prod(n)/k_prod(n))^inv(eta);

    if zbar_prod(n)>1
        l_prod(n)=beta/psi*1/c_prod(n)*(1-al)/(1)*Theta*x_prod(n)^(al/eta)*k_prod(n)^(al*(eta-1)/eta);
        l_prod(n)=l_prod(n)^(inv(al+phi));

        y_prod(n)=Theta*A*x_prod(n)^(al/eta)*k_prod(n)^(al*(eta-1)/eta)*l_prod(n)^(1-al);

    end

    if zbar_prod(n)<=1
        zbar_prod(n)=1;

        l_prod(n)=(1-al)*inv(psi*c_prod(n))*A*((eta)/(eta-1)*lambda*x_prod(n))^al;
        l_prod(n)=l_prod(n)^(inv(al+phi));

        y_prod(n)=A*l_prod(n)^(1-al)*((eta)/(eta-1)*lambda*x_prod(n))^al;

    end


    %     l_prod(n)=beta/psi*1/c_prod(n)*(1-al)/(1)*Theta*x_prod(n)^(al/eta)*k_prod(n)^(al*(eta-1)/eta);
    %     l_prod(n)=l_prod(n)^(inv(al+phi));
    %
    %     y_prod(n)=Theta*x_prod(n)^(al/eta)*k_prod(n)^(al*(eta-1)/eta)*l_prod(n)^(1-al);
    w_prod(n)=(1-al)/(1)*y_prod(n)/l_prod(n);
    r_prod(n)=al*(eta-1)/((1)*eta)*y_prod(n)/k_prod(n);
    %     g_prod(n)=w_prod(n)*l_prod(n)+r_prod(n)*k_prod(n);
    bdot=w_prod(n)*l_prod(n)+r_prod(n)*b_prod(n)-c_prod(n);
    xdot=al/eta*y_prod(n)+(r_prod(n)-delta)*x_prod(n);
    b_prod(n+1)=b_prod(n)+bdot*dt;
    x_prod(n+1)=x_prod(n)+xdot*dt;
    k_prod(n+1)=x_prod(n+1)+b_prod(n+1);

end

save PROD_dynamics.mat

Points=50;

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
 

% x_prod=x_prod(discretizzazione);
% b_prod=b_prod(discretizzazione);
% k_prod=k_prod(discretizzazione);
% y_prod=y_prod(discretizzazione);
% r_prod=r_prod(discretizzazione);
% l_prod=l_prod(discretizzazione);
% w_prod=w_prod(discretizzazione);
% zbar_prod=zbar_prod(discretizzazione);
% c_prod=c_prod(discretizzazione);

%% Plots

% set(gca, 'Units', 'centimeters', 'Position', [2, 2, 25, 10]);

ylim padded

close all

set(gca,'Color','none')
g1 = figure('Name','k_prod');
% plot(discretizzazione,k_prod,'b--','linewidth' , 1)
plot(k_prod,'b--','linewidth' , 1)
% hold on
% plot(CUBueraIR,'k-o','linewidth' , 1)
%h = yline(0, 'k', 'LineWidth', 1, 'HandleVisibility','off');
hold on
plot(discretizzazione,k,'k-o','linewidth' , 1)
hold on
% plot(discretizzazione,k_low,'-*','Color',[251, 86, 7]/256,'linewidth' , 1)

legend('Ambiguity','IM','Interpreter','latex','Location','best')

%xlabel('$t$','Interpreter','latex')
%%title('Producing entrepreneur consumption, IRF','Interpreter','latex')
xlim([0 max(discretizzazione)])
% ylim([1.1*min(k_prod-k)-.1 1.1*max(k_prod-k)+.1])

Ymax=abs(max(max(k_prod),max(k)));
Ymin=abs(min(min(k_prod),min(k)));
YL=Ymin+Ymax;
YLpercentage=0.01*YL;

%ylim([Ymin-YLpercentage Ymax+YLpercentage])

set(gca,'Color','none')
% set(gca,'XTick',[])
ylim padded

print(g1,'k_prod','-depsc',  '-painters','-r600')

close all

set(gca,'Color','none')
g1 = figure('Name','b_prod');
% plot(discretizzazione,b_prod,'b--','linewidth' , 1)
plot(b_prod,'b--','linewidth' , 1)
hold on
plot(discretizzazione,b,'k-o','linewidth' , 1)
hold on
% plot(discretizzazione,b_low,'-*','Color',[251, 86, 7]/256,'linewidth' , 1)
%h = yline(0, 'k', 'LineWidth', 1, 'HandleVisibility','off');
% legend('$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=A$','$\{ A_t \}_{t > 0}=\underline{A}$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','Interpreter','latex','Location','best')
% legend('$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=A$','Interpreter','latex','Location','best')

legend('Ambiguity','IM','Interpreter','latex','Location','best')

%xlabel('$t$','Interpreter','latex')
%%title('Producing entrepreneur consumption, IRF','Interpreter','latex')
xlim([0 max(discretizzazione)])
% ylim([1.1*min(b_prod-b)-.1 1.1*max(b_prod-b)+.1])

Ymax=abs(max(max(b_prod),max(b)));
Ymin=abs(min(min(b_prod),min(b)));
YL=Ymin+Ymax;
YLpercentage=0.01*YL;

%ylim([Ymin-YLpercentage Ymax+YLpercentage])


ylim padded

set(gca,'Color','none')
% set(gca,'XTick',[])

print(g1,'b_prod','-depsc',  '-painters','-r600')

close all


set(gca,'Color','none')
g1 = figure('Name','x_prod');
% plot(discretizzazione,x_prod,'b--','linewidth' , 1)
plot(x_prod,'b--','linewidth' , 1)
hold on
plot(discretizzazione,x,'k-o','linewidth' , 1)

Ymax=abs(max(max(x_prod),max(x)));
Ymin=abs(min(min(x_prod),min(x)));
YL=Ymin+Ymax;
YLpercentage=0.01*YL;

%ylim([Ymin-YLpercentage Ymax+YLpercentage])

% plot(discretizzazione,x_low,'-*','Color',[251, 86, 7]/256,'linewidth' , 1)
%h = yline(0, 'k', 'LineWidth', 1, 'HandleVisibility','off');

% legend('$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=A$','$\{ A_t \}_{t > 0}=\underline{A}$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','Interpreter','latex','Location','best')

% legend('$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=A$','Interpreter','latex','Location','best')

legend('Ambiguity','IM','Interpreter','latex','Location','best')

%xlabel('$t$','Interpreter','latex')
%%title('Producing entrepreneur consumption, IRF','Interpreter','latex')
xlim([0 max(discretizzazione)])
% ylim([1.1*min(x_prod-x)-.1 1.1*max(x_prod-x)+.1])
set(gca,'Color','none')
% set(gca,'XTick',[])

%xlabel('$t$','Interpreter','latex')
%%title('Producing entrepreneur consumption, IRF','Interpreter','latex')
xlim([0 max(discretizzazione)])
% ylim([1.1*min(b_prod-b)-.1 1.1*max(b_prod-b)+.1])
set(gca,'Color','none')
% set(gca,'XTick',[])
ylim padded

print(g1,'x_prod','-depsc',  '-painters','-r600')

close all

set(gca,'Color','none')
g1 = figure('Name','c_prod');
% plot(discretizzazione,c_prod,'b--','linewidth' , 1)
plot(c_prod,'b--','linewidth' , 1)
hold on
plot(discretizzazione,c,'k-o','linewidth' , 1)
hold on
% plot(discretizzazione,c_low,'-*','Color',[251, 86, 7]/256,'linewidth' , 1)
%h = yline(0, 'k', 'LineWidth', 1, 'HandleVisibility','off');

% legend('$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=A$','$\{ A_t \}_{t > 0}=\underline{A}$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','Interpreter','latex','Location','best')
% legend('$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=A$','Interpreter','latex','Location','best')

legend('Ambiguity','IM','Interpreter','latex','Location','best')

%xlabel('$t$','Interpreter','latex')
%%title('Producing entrepreneur consumption, IRF','Interpreter','latex')
xlim([0 max(discretizzazione)])
% ylim([1.1*min(c_prod-c)-.1 1.1*max(c_prod-c)+.1])

Ymax=abs(max(max(c_prod),max(c)));
Ymin=abs(min(min(c_prod),min(c)));
YL=Ymin+Ymax;
YLpercentage=0.01*YL;

%ylim([Ymin-YLpercentage Ymax+YLpercentage])



set(gca,'Color','none')
% set(gca,'XTick',[])
ylim padded

print(g1,'c_prod','-depsc',  '-painters','-r600')

close all

set(gca,'Color','none')
g1 = figure('Name','r_prod');
% plot(discretizzazione,r_prod,'b--','linewidth' , 1)
plot(r_prod,'b--','linewidth' , 1)
hold on
plot(discretizzazione,r,'k-o','linewidth' , 1)
hold on
% plot(discretizzazione,r_low,'-*','Color',[251, 86, 7]/256,'linewidth' , 1)
%h = yline(0, 'k', 'LineWidth', 1, 'HandleVisibility','off');

% legend('$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=A$','$\{ A_t \}_{t > 0}=\underline{A}$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','Interpreter','latex','Location','best')
% legend('$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=A$','Interpreter','latex','Location','best')

legend('Ambiguity','IM','Interpreter','latex','Location','best')

%xlabel('$t$','Interpreter','latex')
%%title('Producing entrepreneur consumption, IRF','Interpreter','latex')
xlim([0 max(discretizzazione)])
% ylim([1.1*min(r_prod-r)-.1 1.1*max(r_prod-r)+.1])

Ymax=abs(max(max(r_prod),max(r)));
Ymin=abs(min(min(r_prod),min(r)));
YL=Ymin+Ymax;
YLpercentage=0.01*YL;

%ylim([Ymin-YLpercentage Ymax+YLpercentage])



set(gca,'Color','none')
% set(gca,'XTick',[])
ylim padded

% 
% axes('position',[.6 .52 .22 .22])  % [left bottom width height]
% %
% % box on % put box around new pair of axes
% %
% plot(discretizzazione(1:min(find(r>r_prod))+2),r_prod(1:min(find(r>r_prod))+2),'b--','linewidth' , 1) % plot on new axes
% hold on
% plot(discretizzazione(1:min(find(r>r_prod))+2),r(1:min(find(r>r_prod))+2),'k-o','linewidth' , 1)     % plot on new axes
% hold on
% plot(r_low(1:min(find(r>r_prod))+1000),'-*','Color',[251, 86, 7]/256,'linewidth' , 1)     % plot on new axes
set(gca,'Color','none')

print(g1,'r_prod','-depsc',  '-painters','-r600')


close all

set(gca,'Color','none')
g1 = figure('Name','zbar_prod');
% plot(discretizzazione,zbar_prod,'b--','linewidth' , 1)
plot(zbar_prod,'b--','linewidth' , 1)
hold on
plot(discretizzazione,zbar,'k-o','linewidth' , 1)
hold on
% plot(discretizzazione,zbar_low,'-*','Color',[251, 86, 7]/256,'linewidth' , 1)
%h = yline(0, 'k', 'LineWidth', 1, 'HandleVisibility','off');

% legend('$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=A$','$\{ A_t \}_{t > 0}=\underline{A}$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','Interpreter','latex','Location','best')
% legend('$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=A$','Interpreter','latex','Location','best')

legend('Ambiguity','IM','Interpreter','latex','Location','best')

%xlabel('$t$','Interpreter','latex')
%%title('Producing entrepreneur consumption, IRF','Interpreter','latex')
xlim([0 max(discretizzazione)])
% ylim([1.1*min(r_prod-r)-.1 1.1*max(r_prod-r)+.1])

Ymax=abs(max(max(zbar_prod),max(zbar)));
Ymin=abs(min(min(zbar_prod),min(zbar)));
YL=Ymin+Ymax;
YLpercentage=0.01*YL;

%ylim([Ymin-YLpercentage Ymax+YLpercentage])



set(gca,'Color','none')
% set(gca,'XTick',[])
ylim padded

print(g1,'zbar_prod','-depsc',  '-painters','-r600')

close all

set(gca,'Color','none')
g1 = figure('Name','w_prod');
% plot(discretizzazione,w_prod,'b--','linewidth' , 1)
plot(w_prod,'b--','linewidth' , 1)
hold on
plot(discretizzazione,w,'k-o','linewidth' , 1)
hold on
% plot(discretizzazione,w_low,'-*','Color',[251, 86, 7]/256,'linewidth' , 1)
%h = yline(0, 'k', 'LineWidth', 1, 'HandleVisibility','off');

% legend('$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=A$','$\{ A_t \}_{t > 0}=\underline{A}$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','Interpreter','latex','Location','best')

legend('$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=A$','Interpreter','latex','Location','best')

legend('Ambiguity','IM','Interpreter','latex','Location','best')

%xlabel('$t$','Interpreter','latex')
%%title('Producing entrepreneur consumption, IRF','Interpreter','latex')
xlim([0 max(discretizzazione)])
% ylim([1.1*min(r_prod-r)-.1 1.1*max(r_prod-r)+.1])
set(gca,'Color','none')
% set(gca,'XTick',[])

% axes('position',[.6 .20 .22 .22])  % [left bottom width height]
% 
% % box on % put box around new pair of axes
% 
% plot(discretizzazione(1:min(find(w_prod>w))+2),w_prod(1:min(find(w_prod>w))+2),'b--','linewidth' , 1) % plot on new axes
% hold on
% plot(discretizzazione(1:min(find(w_prod>w))+2),w(1:min(find(w_prod>w))+2),'k-o','linewidth' , 1)     % plot on new axes
% hold on
% plot(w_low(1:min(find(w_prod>w))+1000),'-*','Color',[251, 86, 7]/256,'linewidth' , 1)     % plot on new axes
set(gca,'Color','none')

Ymax=abs(max(max(w_prod),max(w)));
Ymin=abs(min(min(w_prod),min(w)));
YL=Ymin+Ymax;
YLpercentage=0.01*YL;

%ylim([Ymin-YLpercentage Ymax+YLpercentage])


ylim padded

print(g1,'w_prod','-depsc',  '-painters','-r600')

close all

set(gca,'Color','none')
g1 = figure('Name','y_prod');
% plot(discretizzazione,y_prod,'b--','linewidth' , 1)
plot(y_prod,'b--','linewidth' , 1)
hold on
plot(discretizzazione,y,'k-o','linewidth' , 1)
hold on
% plot(discretizzazione,y_low,'-*','Color',[251, 86, 7]/256,'linewidth' , 1)
% %h = yline(0, 'k', 'LineWidth', 1, 'HandleVisibility','off');

% legend('$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=A$','$\{ A_t \}_{t > 0}=\underline{A}$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','Interpreter','latex','Location','best')
% legend('$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=A$','Interpreter','latex','Location','best')

legend('Ambiguity','IM','Interpreter','latex','Location','best')

%xlabel('$t$','Interpreter','latex')
%%title('Producing entrepreneur consumption, IRF','Interpreter','latex')
xlim([0 max(discretizzazione)])

Ymax=abs(max(max(y_prod),max(y)));
Ymin=abs(min(min(y_prod),min(y)));
YL=Ymin+Ymax;
YLpercentage=0.01*YL;

%ylim([Ymin-YLpercentage Ymax+YLpercentage])

set(gca,'Color','none')
% set(gca,'XTick',[])
ylim padded





print(g1,'y_prod','-depsc',  '-painters','-r600')

close all

set(gca,'Color','none')
g1 = figure('Name','l_prod');
% plot(discretizzazione,l_prod,'b--','linewidth' , 1)
plot(l_prod,'b--','linewidth' , 1)
hold on
plot(discretizzazione,l,'k-o','linewidth' , 1)
hold on

% plot(discretizzazione,l_low,'-*','Color',[251, 86, 7]/256,'linewidth' , 1)

% %h = yline(0, 'k', 'LineWidth', 1, 'HandleVisibility','off');

% legend('$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=A$','$\{ A_t \}_{t > 0}=\underline{A}$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','Interpreter','latex','Location','best')

legend('$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=\underline{A}$','$\{ A_t \}_{t \geq 0}=A$, $\{ \underline{A}_t \}_{t > 0}=A$','Interpreter','latex','Location','best')

legend('Ambiguity','IM','Interpreter','latex','Location','best')


%xlabel('$t$','Interpreter','latex')
%%title('Producing entrepreneur consumption, IRF','Interpreter','latex')
xlim([0 max(discretizzazione)])

Ymax=abs(max(max(l_prod),max(l)));
Ymin=abs(min(min(l_prod),min(l)));
YL=Ymin+Ymax;
YLpercentage=0.01*YL;

%ylim([Ymin-YLpercentage Ymax+YLpercentage])


set(gca,'Color','none')
% set(gca,'XTick',[])
ylim padded

print(g1,'l_prod','-depsc',  '-painters','-r600')

