rep=100;
Nday=63;
S0=100;K=105;sigma=0.2;r=0.01;T=1;
sigmag=0.1;
dt = 1/251;
S=S0*ones(Nday+1,rep);
for i=1:Nday
    rng(i);
S(i+1,:)=S(i,:).*exp((r-0.5*sigmag^2)*dt+sigmag*sqrt(dt)*normrnd(0,1,1,rep));
end
portfolio_ESK=zeros(rep,Nday-1);
portfolio_SK=zeros(rep,Nday-1);
portfolio_true=zeros(rep,Nday-1);
SS=zeros(Nday+1,rep);
for i=1:rep
    rng(i+rep);
[p_ESK,p_SK,p_true,S1]=derivative_gap(Nday,S(:,i));
portfolio_ESK(i,:)=p_ESK;
portfolio_SK(i,:)=p_SK;
portfolio_true(i,:)=p_true;
SS(:,i)=S1;
end

var_port_ESK = sqrt(var(portfolio_ESK));
var_port_SK = sqrt(var(portfolio_SK));
var_port_true = sqrt(var(portfolio_true));
plot(1:Nday-1,var_port_ESK,'r*-',1:Nday-1,var_port_SK,'bo--',1:Nday-1,var_port_true,'ks:','LineWidth',2)
xlabel('Trading days')
ylabel('Standard deviation of Profit & Loss')
set(gca,'FontSize',16);
legend('GESK','SK','True')
save('63days_var')
