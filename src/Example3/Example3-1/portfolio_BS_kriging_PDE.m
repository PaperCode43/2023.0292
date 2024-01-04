%% consider the portfolio has one Basket option
function [RMSE_V,RMSE_theta,ARMSE_delta,ARMSE_gamma]=portfolio_BS_kriging_PDE(X,poly_d,seed,M,M0)
%tic
r=0.05;%interest rate
sigma=0.2;%volatility
T=1;%maturity
S0=101:140;%initial values
K=81:2:159;
nl=length(S0);

N=5000;% number of observations to price options
Nt=500;% number of testing points
Nd=1000;% number of dicretized points for appoximating the integral and expectation in the pdf-objective function
pd=sobolset(2);
%% specify the random number

Y_true=zeros(Nt,nl);delta_true=zeros(Nt,nl);gamma_true=zeros(Nt,nl);theta_true=zeros(Nt,nl);
f_ESK=zeros(Nt,nl);f_ESK_pde=zeros(Nt,nl);
delta_ESK=zeros(Nt,nl);delta_ESK_pde=zeros(Nt,nl);
gamma_ESK=zeros(Nt,nl);gamma_ESK_pde=zeros(Nt,nl);
theta_ESK=zeros(Nt,nl);theta_ESK_pde=zeros(Nt,nl);
parfor i=1:nl
    rng(100*seed+i);
    S_taumax=S0(i)*1.2;t_max=T; %range of design points [5,25]*[0,T]
    S_taumin=S0(i)*0.8;t_min=0;
    S_taumind=S0(i)*0.8;S_taumaxd=S0(i)*1.2;
    t_mind=0;t_maxd=T;
%% generate testing points 
    u=pd(2+(i-1)*Nt:i*Nt+1,:);
    Stest = S_taumin+(S_taumax-S_taumin).*u(:,1);
    ttest=t_min+(t_max-t_min).*u(:,2);
    Xtest=[Stest,ttest];
    [Y1_true,delta1_true,gamma1_true,theta1_true]=BS_price(Xtest(:,1),r,sigma,T-Xtest(:,2),K(i));
    Y_true(:,i)=Y1_true;delta_true(:,i)=delta1_true;gamma_true(:,i)=gamma1_true;theta_true(:,i)=theta1_true;
    %% Generate design points
    p=lhsdesign(M,2);
    p2=lhsdesign(M0,1);
    S_b=(S_taumind+(S_taumaxd-S_taumind)*p2);
    S0d=[(S_taumind+(S_taumaxd-S_taumind)*p(:,1));S_b];
    t0=[t_mind+(t_maxd-t_mind)*p(:,2);ones(M0,1)*T];
    [Y1,vY1,dS1,Vd1,Vmatrix1]=BS_MC(S0d,r,sigma,T-t0,N,K(i));
    X=[S0d,t0];
    [Wm_price1,Wm_delta1,Wm_gamma1,Wm_theta1,Wm_beta1,~]...
         =BS_portfolio_kriging_PDE_grad(r,sigma,S0(i),T,Nd,X,Y1,vY1,Vmatrix1,Xtest,poly_d,dS1);

    B1=polybasis(X,poly_d);
    Btest=polybasis(Xtest,poly_d);
    d_Btest_only=polybasis(Xtest,poly_d);
    d_Btest_delta=d_polybasis(Xtest,poly_d,2);
    d_Btest_gamma=d_Btest_delta;% all zero
    d_Btest_theta=d_Btest_delta;% all zero
    d_B_dS=d_polybasis(X,poly_d,1);

    betahat1=Wm_beta1.SK;
    beta_pde1=Wm_beta1.SKpde;
    dfB_dS1=d_B_dS*betahat1;
    dYnew_dS1=dS1 - dfB_dS1;
    dfB_dSpde1=d_B_dS*beta_pde1;
    dYnew_dSpde1=dS1 - dfB_dSpde1;
    %% SK surface 
    Ynew1=Y1 - B1*betahat1;
    f_SK1=Wm_price1.SK'*Ynew1+Btest*betahat1;
    %% SKpde surface
    Ynew_pde1=Y1 - B1*beta_pde1;
    f_SK_pde1=Wm_price1.SKpde'*Ynew_pde1+Btest*Wm_beta1.SK;
    %% ESK surface
    f_ESK(:,i)=Wm_price1.ESK'*[Ynew1;dYnew_dS1]+Btest*Wm_beta1.ESK;
    %% ESKpde surface
    f_ESK_pde(:,i)=Wm_price1.ESKpde'*[Ynew_pde1;dYnew_dSpde1]+Btest*Wm_beta1.ESKpde;
    %% Greeks
    %% SK
    delta_SK1 = Wm_delta1.SK'*Ynew1 + d_Btest_delta*betahat1;
    gamma_SK1 = Wm_gamma1.SK'*Ynew1 + d_Btest_gamma*betahat1;    
    theta_SK1 = Wm_theta1.SK'*Ynew1 + d_Btest_theta*betahat1;
    %% SKpde
    delta_SK_pde1 = Wm_delta1.SKpde'*Ynew_pde1 + d_Btest_delta*Wm_beta1.SKpde;
    gamma_SK_pde1 = Wm_gamma1.SKpde'*Ynew_pde1 + d_Btest_gamma*Wm_beta1.SKpde;
    theta_SK_pde1 = Wm_theta1.SKpde'*Ynew_pde1 + d_Btest_theta*Wm_beta1.SKpde;
    %% ESK
    delta_ESK(:,i) = Wm_delta1.ESK'*[Ynew1;dYnew_dS1] + d_Btest_delta*betahat1;
    gamma_ESK(:,i) = Wm_gamma1.ESK'*[Ynew1;dYnew_dS1] + d_Btest_gamma*betahat1;
    theta_ESK(:,i) = Wm_theta1.ESK'*[Ynew1;dYnew_dS1] + d_Btest_theta*betahat1;
    %% ESKpde
    delta_ESK_pde(:,i) = Wm_delta1.ESKpde'*[Ynew_pde1;dYnew_dSpde1] + d_Btest_delta*Wm_beta1.ESKpde;
    gamma_ESK_pde(:,i) = Wm_gamma1.ESKpde'*[Ynew_pde1;dYnew_dSpde1] + d_Btest_gamma*Wm_beta1.ESKpde;
    theta_ESK_pde(:,i) = Wm_theta1.ESKpde'*[Ynew_pde1;dYnew_dSpde1] + d_Btest_theta*Wm_beta1.ESKpde;
end
Port_true=sum(Y_true,2);Port_T_true=sum(theta_true,2);
Port_ESK = sum(f_ESK,2);Port_ESK_pde=sum(f_ESK_pde,2);
Port_T_ESK = sum(theta_ESK,2);Port_T_ESK_pde=sum(theta_ESK_pde,2);

RMSE_V.ESK = sqrt(mean((Port_ESK - Port_true).^2));
RMSE_V.ESKpde = sqrt(mean((Port_ESK_pde - Port_true).^2));
RMSE_theta.ESK = sqrt(mean((Port_T_ESK - Port_T_true).^2));
RMSE_theta.ESKpde = sqrt(mean((Port_T_ESK_pde - Port_T_true).^2));


RMSE_delta.ESK=zeros(nl,1);RMSE_delta.ESKpde=zeros(nl,1);
RMSE_gamma.ESK=zeros(nl,1);RMSE_gamma.ESKpde=zeros(nl,1);
for i=1:nl
    RMSE_delta.ESK(i) = sqrt(mean((delta_ESK(:,i) - delta_true(:,i)).^2));
    RMSE_delta.ESKpde(i) = sqrt(mean((delta_ESK_pde(:,i) - delta_true(:,i)).^2));
    RMSE_gamma.ESK(i) = sqrt(mean((gamma_ESK(:,i) - gamma_true(:,i)).^2));
    RMSE_gamma.ESKpde(i) = sqrt(mean((gamma_ESK_pde(:,i) - gamma_true(:,i)).^2));
end

ARMSE_delta.ESK=mean(RMSE_delta.ESK);
ARMSE_delta.ESKpde=mean(RMSE_delta.ESKpde);
ARMSE_gamma.ESK = mean(RMSE_gamma.ESK);
ARMSE_gamma.ESKpde = mean(RMSE_gamma.ESKpde);









function basis=polybasis(X,deg)
basis=ones(size(X,1),1);
for i=1:deg
    basis=[basis,X.^i];
end

function d_basis=d_polybasis(X,deg,partial_idx)
if deg==0
    d_basis=zeros(size(X,1),1);
elseif deg==1
    Xn=zeros(size(X));
    Xn(:,partial_idx) = ones(size(X,1),1);
    d_basis=[zeros(size(X,1),1),Xn];
else
Xn1=zeros(size(X));
Xn1(:,partial_idx) = ones(size(X,1),1);
d_basis=[zeros(size(X,1),1),Xn1];
Xn=zeros(size(X));
Xn(:,partial_idx) = X(:,partial_idx);
for i=2:deg
    d_basis=[d_basis,(i)*Xn.^(i-1)];
end
end
    
    

