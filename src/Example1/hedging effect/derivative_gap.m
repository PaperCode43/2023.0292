function [gap_ESK,gap_SK,gap_true,S] = derivative_gap(Nday,S)
S0=100;K=105;sigma=0.2;r=0.01;T=1;M=20000;t=0;poly_d=0;Tru_MC_idx=1;
dt = 1/251;

%% generate Ndays stock price
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
X=S(i)*[0.96,0.99,1.03]';
Xtest=S(i+1);
[f_true(i),df_true(i),f_ESK(i),df_ESK(i),f_SK(i),df_SK(i)]...%% (i) is time series
    =KrigEuro_cokriging_delta_1(X,r,sigma,T-i*dt,K,M,t,poly_d,Tru_MC_idx,Xtest);
%[Y_true,dY_true,f_enhanced,df_enhanced,f_SK,df_SK]=KrigEuro_cokriging_delta(X,r,sigma,T-dt*i,K,M,t,poly_d,Tru_MC_idx);
end
option_change_ESK=zeros(Nday,1);option_change_SK=zeros(Nday,1);option_change_true=zeros(Nday,1);
delta_hedge_ESK=zeros(Nday,1);delta_hedge_SK=zeros(Nday,1);delta_hedge_true=zeros(Nday,1);
%portvalue_SK=zeros(Nday,1);portvalue_ESK=zeros(Nday,1);portvalue_true=zeros(Nday,1);
for i=2:Nday
    option_change_ESK(i) = f_ESK(i) - f_ESK(i-1);
    option_change_SK(i) = f_SK(i) - f_SK(i-1);
    option_change_true(i) = f_true(i) - f_true(i-1);
    delta_hedge_ESK(i) = (S(i+1) - S(i)) * df_ESK(i); % for S(1) corresponse S0, S(2) =>S1 => f_SK(1)
    delta_hedge_SK(i) = (S(i+1) - S(i)) * df_SK(i);
    delta_hedge_true(i) = (S(i+1) - S(i)) * df_true(i);
    %portvalue_SK(i) = -f_SK(i) + df_SK(i)*S(i+1);
    %portvalue_ESK(i) = -f_ESK(i) + df_ESK(i)*S(i+1);
    %portvalue_true(i) = -f_true(i) + df_true(i)*S(i+1);
end
gap_ESK = (option_change_ESK(2:end) - delta_hedge_ESK(2:end))*1000;
gap_SK = (option_change_SK(2:end) - delta_hedge_SK(2:end))*1000;
gap_true = (option_change_true(2:end) - delta_hedge_true(2:end))*1000;
plot(1:Nday-1,gap_ESK,'r*-',1:Nday-1,gap_SK,'bo--',1:Nday-1,gap_true,'ks:','LineWidth',2)
xlabel('Trading days')
ylabel('Profit & Loss')
set(gca,'FontSize',16);
legend('GESK','SK','True')
%plot(1:Nday-1,portvalue_ESK(2:end),'r*-',1:Nday-1,portvalue_SK(2:end),'bo--',1:Nday-1,portvalue_true(2:end),'ks--','LineWidth',2)
