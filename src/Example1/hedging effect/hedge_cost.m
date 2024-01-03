rep=100;
Nday=63;
S0=100;K=105;sigma=0.2;r=0.01;T=1;
dt = 1/251;
sigmag=0.1;

S=S0*ones(Nday+1,rep);
S_adj=S0*ones(Nday+1,rep);
for i=1:Nday
    rng(i);
S(i+1,:)=S(i,:).*exp((r-0.5*sigmag^2)*dt+sigmag*sqrt(dt)*normrnd(0,1,1,rep));
S_adj(i+1,:)=S_adj(i,:).*exp((r-0.5*sigmag^2)*dt+sigmag*sqrt(dt)*normrnd(0,1,1,rep));
end

cost_ESK=zeros(rep,1);cost_SK=zeros(rep,1);cost_SK_SE=zeros(rep,1);
parfor i=1:rep
    rng(i+rep);
[cost_ESK(i),cost_SK(i),cost_SK_SE(i),~] = hedge_cost_once(Nday,S(:,i),S_adj(:,i));
end

figure(2)
hold on
boxplot([cost_SK_SE,cost_SK,cost_ESK],'Labels',{'SK same effect','SK same frequency','GESK'})
ylabel('Hedging cost')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off
save('63_period_hedging_cost_fitconst_10design')
