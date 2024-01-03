function [cost_ESK,cost_SK,cost_SK_SE,S] = hedge_cost_once(Nday,S,S_adj)
S0=100;K=105;sigma=0.2;r=0.01;T=1;M=5000;t=0;poly_d=0;Tru_MC_idx=0;%final setting
dt = 1/251;
%sigmag=0.4;
% %% generate Ndays stock price
% S=S0*ones(Nday+1,1);
% for i=1:Nday
% S(i+1)=S(i)*exp((r-0.5*sigmag^2)*dt+sigmag*sqrt(dt)*normrnd(0,1,1,1));
% end
%% consider Ndays hedging, generate ESK option price and ESK delta,  
%%% and SK option price and SK delta on these days
f_ESK = zeros(Nday,1);
df_ESK = zeros(Nday,1);
f_SK = zeros(Nday,1);
df_SK = zeros(Nday,1);
f_true=zeros(Nday,1);
df_true=zeros(Nday,1);
for i=1:Nday
    u=0.8+0.4*lhsdesign(10,1);%final setting
X=S(i)*u;
%X=S(i)*[0.96,0.99,1.03]';
Xtest=S(i+1);
[f_true(i),df_true(i),f_ESK(i),df_ESK(i),f_SK(i),df_SK(i)]...%% (i) is time series
    =KrigEuro_cokriging_delta_1(X,r,sigma,T-i*dt,K,M,t,poly_d,Tru_MC_idx,Xtest);
%[Y_true,dY_true,f_enhanced,df_enhanced,f_SK,df_SK]=KrigEuro_cokriging_delta(X,r,sigma,T-dt*i,K,M,t,poly_d,Tru_MC_idx);
end
adjust_SK=zeros(Nday,1);
adjust_extra_SK=zeros(Nday,1);
adjust_ESK=zeros(Nday,1);
for i=2:Nday
    adjust_ESK(i) = abs(df_ESK(i) - df_ESK(i-1))*S(i+1)*1000*0.001;
    adjust_SK(i) = abs(df_SK(i) - df_SK(i-1))*S(i+1)*1000*0.001;
    adjust_extra_SK(i) = abs(df_ESK(i-1)-df_SK(i-1))*S_adj(i+1)*1000*0.001;
    %portvalue_SK(i) = -f_SK(i) + df_SK(i)*S(i+1);
    %portvalue_ESK(i) = -f_ESK(i) + df_ESK(i)*S(i+1);
    %portvalue_true(i) = -f_true(i) + df_true(i)*S(i+1);
end
cost_ESK = sum(adjust_ESK);
cost_SK_SE = sum(adjust_SK)+sum(adjust_extra_SK);
cost_SK = sum(adjust_SK);
% plot(1:Nday-1,gap_ESK,'r*-',1:Nday-1,gap_SK,'bo--',1:Nday-1,gap_true,'ks:','LineWidth',2)
% xlabel('Trading days')
% ylabel('Profit & Loss')
% set(gca,'FontSize',16);
% legend('GESK','SK','True')
%plot(1:Nday-1,portvalue_ESK(2:end),'r*-',1:Nday-1,portvalue_SK(2:end),'bo--',1:Nday-1,portvalue_true(2:end),'ks--','LineWidth',2)
