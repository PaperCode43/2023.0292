function [portfolio_ESK,portfolio_SK,S]=Portfolio_gap(Nday)
sigma_v=[0.2636 0.2625 0.4012 0.2842 0.4660];
nu_v = [0.0387 0.0355 0.0394 0.0017 0.0933];
theta_v = [-0.5185 -0.8288 -1.2344 -2.6984 -1.1459];
S0_v=[204.47 144.96 113.47 144.85 120.51];
r_n=0.0273;
K_v = [204,140,110,150,125];% one at-the-money, two in the money, two out of money 
T=1;
n=12;
M=40000;
S=ones(Nday+1,1)*S0_v;
sigmag=0.1;
dt=T/(30*n);
ZZ=normrnd(0,1,Nday,5);

for i=1:Nday
S(i+1,:)=S(i,:).*exp((r_n-0.5*sigmag^2)*dt+sigmag*sqrt(dt)*ZZ(i,:));
end
%%% hedge gap
[AAPL_gap_Asian_ESK,AAPL_gap_Asian_SK,AAPL_gap_LB_ESK,AAPL_gap_LB_SK] = hedge_gap(sigma_v(1),nu_v(1),theta_v(1),S0_v(1),r_n,K_v(1),T,n,Nday,M,S(:,1));
[FB_gap_Asian_ESK,FB_gap_Asian_SK,FB_gap_LB_ESK,FB_gap_LB_SK] = hedge_gap(sigma_v(2),nu_v(2),theta_v(2),S0_v(2),r_n,K_v(2),T,n,Nday,M,S(:,2));
[NFL_gap_Asian_ESK,NFL_gap_Asian_SK,NFL_gap_LB_ESK,NFL_gap_LB_SK] = hedge_gap(sigma_v(3),nu_v(3),theta_v(3),S0_v(3),r_n,K_v(3),T,n,Nday,M,S(:,3));
[ALI_gap_Asian_ESK,ALI_gap_Asian_SK,ALI_gap_LB_ESK,ALI_gap_LB_SK] = hedge_gap(sigma_v(4),nu_v(4),theta_v(4),S0_v(4),r_n,K_v(4),T,n,Nday,M,S(:,4));
[TSL_gap_Asian_ESK,TSL_gap_Asian_SK,TSL_gap_LB_ESK,TSL_gap_LB_SK] = hedge_gap(sigma_v(5),nu_v(5),theta_v(5),S0_v(5),r_n,K_v(5),T,n,Nday,M,S(:,5));

portfolio_ESK=AAPL_gap_Asian_ESK+AAPL_gap_LB_ESK+FB_gap_Asian_ESK+FB_gap_LB_ESK+NFL_gap_Asian_ESK+NFL_gap_LB_ESK...
               +ALI_gap_Asian_ESK+ALI_gap_LB_ESK+TSL_gap_Asian_ESK+TSL_gap_LB_ESK;
portfolio_SK=AAPL_gap_Asian_SK+AAPL_gap_LB_SK+FB_gap_Asian_SK+FB_gap_LB_SK+NFL_gap_Asian_SK+NFL_gap_LB_SK...
               +ALI_gap_Asian_SK+ALI_gap_LB_SK+TSL_gap_Asian_SK+TSL_gap_LB_SK;
%plot(1:29,portfolio_ESK,'r*-',1:29,portfolio_SK,'bo--',1:29,gap_true,'kv:')
plot(1:Nday-1,portfolio_ESK,'r*-',1:Nday-1,portfolio_SK,'bo--','LineWidth',2)
%save('63daydata')
