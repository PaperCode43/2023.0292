function [gap_Asian_ESK,gap_Asian_SK,gap_LB_ESK,gap_LB_SK] = hedge_gap(sigma_v,nu_v,theta_v,S0_v,r_n,K_v,T,n,Nday,M,S)
%% generate Ndays stocks price
poly_d=0;
% AAPL

%% begin fitting
% AAPL Asian
f_Asian_ESK = zeros(Nday,1);
df_Asian_ESK = zeros(Nday,1);
f_Asian_SK = zeros(Nday,1);
df_Asian_SK = zeros(Nday,1);
f_LB_ESK = zeros(Nday,1);
df_LB_ESK = zeros(Nday,1);
f_LB_SK = zeros(Nday,1);
df_LB_SK = zeros(Nday,1);
for i=1:Nday
X=S(i)*[0.96,0.99,1.04]';
Xtest=S(i+1);
[V_A,VdS0_A]=KrigEuro_cokriging_delta_Asian1(X,K_v,T,i,n,r_n,sigma_v,nu_v,theta_v,M,poly_d,Xtest);
[V_L,VdS0_L]=KrigEuro_cokriging_delta_Lookback1(X,K_v,T,i,n,r_n,sigma_v,nu_v,theta_v,M,poly_d,Xtest);
% alpha is the range precentage of the S0
f_Asian_ESK(i)=V_A.ESK;df_Asian_ESK(i)=VdS0_A.ESK;f_Asian_SK(i)=V_A.SK;df_Asian_SK(i)=VdS0_A.SK;
f_LB_ESK(i)=V_L.ESK;df_LB_ESK(i)=VdS0_L.ESK;f_LB_SK(i)=V_L.SK;df_LB_SK(i)=VdS0_L.SK;
end
option_change_Asian_ESK=zeros(Nday,1);option_change_Asian_SK=zeros(Nday,1);
delta_hedge_Asian_ESK=zeros(Nday,1);delta_hedge_Asian_SK=zeros(Nday,1);
option_change_LB_ESK=zeros(Nday,1);option_change_LB_SK=zeros(Nday,1);
delta_hedge_LB_ESK=zeros(Nday,1);delta_hedge_LB_SK=zeros(Nday,1);
for i=2:Nday
    option_change_Asian_ESK(i) = f_Asian_ESK(i) - f_Asian_ESK(i-1);
    option_change_Asian_SK(i) = f_Asian_SK(i) - f_Asian_SK(i-1);
    delta_hedge_Asian_ESK(i) = (S(i+1) - S(i)) * df_Asian_ESK(i); % for S(1) corresponse S0, S(2) =>S1 => f_SK(1)
    delta_hedge_Asian_SK(i) = (S(i+1) - S(i)) * df_Asian_SK(i);
    option_change_LB_ESK(i) = f_LB_ESK(i) - f_LB_ESK(i-1);
    option_change_LB_SK(i) = f_LB_SK(i) - f_LB_SK(i-1);
    delta_hedge_LB_ESK(i) = (S(i+1) - S(i)) * df_LB_ESK(i); % for S(1) corresponse S0, S(2) =>S1 => f_SK(1)
    delta_hedge_LB_SK(i) = (S(i+1) - S(i)) * df_LB_SK(i);
end
gap_Asian_ESK = option_change_Asian_ESK(2:end) - delta_hedge_Asian_ESK(2:end);
gap_Asian_SK = option_change_Asian_SK(2:end) - delta_hedge_Asian_SK(2:end);
gap_LB_ESK = option_change_LB_ESK(2:end) - delta_hedge_LB_ESK(2:end);
gap_LB_SK = option_change_LB_SK(2:end) - delta_hedge_LB_SK(2:end);


