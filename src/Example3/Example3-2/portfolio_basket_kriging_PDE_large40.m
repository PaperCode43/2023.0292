%% consider the portfolio has one Basket option
function [RMSE_V,RMSE_theta,ARMSE_delta1,ARMSE_delta2,ARMSE_gamma11,ARMSE_gamma12,ARMSE_gamma21,ARMSE_gamma22, RMSE_delta1, RMSE_delta2,RMSE_gamma11,RMSE_gamma12,RMSE_gamma22]=portfolio_basket_kriging_PDE_large40(X,poly_d,seed,M,M0)
%tic
r=0.02;%interest rate
sigma=[0.1*ones(10,1),0.15*ones(10,1);0.2*ones(10,1),0.15*ones(10,1);0.2*ones(10,1),0.2*ones(10,1);0.25*ones(10,1),0.2*ones(10,1)];%volatility
T=1;%maturity
S0=[(10:1:49)',(20:1:59)'];%initial values
K1=8:0.6:31.4;K2=10:0.7:37.3;K3=12:0.8:43.2;%strike prices
K4=16:1.1:58.9;K5=18:1.2:64.8;K6=20:1.3:70.7;
K=[K1',K2',K3',K4',K5',K6'];
Rho=zeros(2,2,40);
for i=1:10
Rho(:,:,i)=[1 0.1;0.1 1];Rho(:,:,10+i)=[1 0.5;0.5 1];Rho(:,:,20+i)=[1 -0.5;-0.5 1];Rho(:,:,30+i)=[1 0; 0 1];
end

N=4000;% number of observations to price options
Nt=500;% number of testing points
Nd=1000;% number of dicretized points for appoximating the integral and expectation in the pdf-objective function
pd=sobolset(3);
%% specify the random number

Y_true=zeros(Nt,40);delta_true=zeros(Nt,2,40);gamma_true=zeros(Nt,4,40);theta_true=zeros(Nt,40);
f_ESK=zeros(Nt,40);f_ESK_pde=zeros(Nt,40);
delta1_ESK=zeros(Nt,40);delta2_ESK=zeros(Nt,40);delta1_ESK_pde=zeros(Nt,40);delta2_ESK_pde=zeros(Nt,40);
gamma11_ESK=zeros(Nt,40);gamma12_ESK=zeros(Nt,40);gamma21_ESK=zeros(Nt,40);gamma22_ESK=zeros(Nt,40);
gamma11_ESK_pde=zeros(Nt,40);gamma12_ESK_pde=zeros(Nt,40);gamma21_ESK_pde=zeros(Nt,40);gamma22_ESK_pde=zeros(Nt,40);
theta_ESK=zeros(Nt,40);theta_ESK_pde=zeros(Nt,40);
Mfinpar=zeros(5,2,40);
parfor i=1:40
    %i
    rng(100*seed+i);
    S_taumax=S0(i,:)*1.5;t_max=T; %range of design points [5,25]*[0,T]
    S_taumin=S0(i,:)*0.5;t_min=0;
    S_taumind=S0(i,:)*0.48;S_taumaxd=S0(i,:)*1.52;
    t_mind=0;t_maxd=1;
    %% generate testing points 
    u=pd(2+(i-1)*Nt:i*Nt+1,:);
    Stest = ones(Nt,1)*S_taumin+(ones(Nt,1)*S_taumax-ones(Nt,1)*S_taumin).*u(:,1:2);
    ttest=t_min+(t_max-t_min).*u(:,3);
    Xtest=[Stest,ttest];
    %[Y_true,delta_true,gamma_true,theta_true]=portfolio_basket_price(Xtest(:,1:2),r,sigma,T-Xtest(:,3),rho,Nt,K1,K2,K3);
    %Y_true is Nt*1, delta_true is Nt*2, gamma_true is Nt*4 (d1d1,d1d2,d2d1,d2d2), theta_true is Nt*1
    [Y1_true,delta1_true,gamma1_true,theta1_true]=portfolio_basket_price2(Xtest(:,1:2),r,sigma(i,:),T-Xtest(:,3),Rho(:,:,i),K1(i),K2(i),K3(i));
    [Y2_true,delta2_true,gamma2_true,theta2_true]=portfolio_basket_price2(Xtest(:,1:2),r,sigma(i,:),T-Xtest(:,3),Rho(:,:,i),K4(i),K5(i),K6(i));
    Y_true(:,i)=10*Y1_true+Y2_true;delta_true(:,:,i)=10*delta1_true+delta2_true;gamma_true(:,:,i)=10*gamma1_true+gamma2_true;theta_true(:,i)=-10*theta1_true-theta2_true;  
%     [Y1_truep,delta1_true,gamma1_true,theta1_true]=portfolio_basket_price2(Xtest(:,1:2),r,sigma(i,:),T-(Xtest(:,3)+0.001),Rho(:,:,i),K1(i),K2(i),K3(i));
%     [Y2_truep,delta2_true,gamma2_true,theta2_true]=portfolio_basket_price2(Xtest(:,1:2),r,sigma(i,:),T-(Xtest(:,3)+0.001),Rho(:,:,i),K4(i),K5(i),K6(i));
%     [Y1_truen,delta1_true,gamma1_true,theta1_true]=portfolio_basket_price2(Xtest(:,1:2),r,sigma(i,:),T-(Xtest(:,3)-0.001),Rho(:,:,i),K1(i),K2(i),K3(i));
%     [Y2_truen,delta2_true,gamma2_true,theta2_true]=portfolio_basket_price2(Xtest(:,1:2),r,sigma(i,:),T-(Xtest(:,3)-0.001),Rho(:,:,i),K4(i),K5(i),K6(i));
%     
%     theta_FD=((10*Y1_truep+Y2_truep)-(10*Y1_truen+Y2_truen))/0.002;  
    %% Generate design points
    p=lhsdesign(M,3);
    p2=lhsdesign(M0,2);
    S_b=(S_taumind+(S_taumaxd-S_taumind).*p2);
    S0d=[(S_taumind+(S_taumaxd-S_taumind).*p(:,1:2));S_b];
    t0=[t_mind+(t_maxd-t_mind)*p(:,3);ones(M0,1)*T];
    [Y,vY,dS1,Vd1,dS2,Vd2,Vmatrix]=basket_MC3(S0d,r,sigma(i,:),T-t0,Rho(:,:,i),N,K1(i),K2(i),K3(i),K4(i),K5(i),K6(i));
    %[Y_trued1p,dS11p,~,~]=portfolio_basket_price2(S0d+[0,0.0000001],r,sigma,T-t0,Rho,K1,K2,K3);
    %[Y_trued2p,dS12p,~,~]=portfolio_basket_price2(S0d+[0,0.0000001],r,sigma,T-t0,Rho,K4,K5,K6);
    %Y_truedp=10*Y_trued1p+Y_trued2p;
    %dS1_truedp=10*dS11p+dS12p;
    % [Y_trued1n,~,~,~]=portfolio_basket_price2(S0-[0,0.0001],r,sigma,T-t0,Rho,K1,K2,K3);
    % [Y_trued2n,~,~,~]=portfolio_basket_price2(S0-[0,0.0001],r,sigma,T-t0,Rho,K4,K5,K6);
    % Y_truedn=10*Y_trued1n+Y_trued2n;
    % delta2FD=(Y_truedp - Y_truedn)/0.0002;
    X=[S0d,t0];
    [Wm_price,Wm_delta1,Wm_delta2,Wm_gamma11,Wm_gamma12,Wm_gamma22,Wm_theta,Wm_beta,finpar]=basket_kriging_PDE_grad(r,sigma(i,:),S0(i,:),K(i,:),T,Rho(:,:,i),Nd,X,Y,vY,Vmatrix,Xtest,i,poly_d,dS1,dS2);
    fp=finpar;
    Mfinpar(:,:,i)=fp;
    betahat=Wm_beta.SK;
    beta_pde=Wm_beta.SKpde;
    B1=polybasis(X,poly_d);
    Btest=polybasis(Xtest,poly_d);
    d_Btest_delta=d_polybasis(Xtest,poly_d,2);
    d_Btest_gamma=d_Btest_delta;% all zero
    d_Btest_theta=d_Btest_delta;% all zero
    d_B_dS1=d_polybasis(X,poly_d,1);
    dfB_dS1=d_B_dS1*betahat;
    dYnew_dS1=dS1 - dfB_dS1;
    dfB_dS1pde=d_B_dS1*beta_pde;
    dYnew_dS1pde=dS1 - dfB_dS1pde;
    d_B_dS2=d_polybasis(X,poly_d,2);
    dfB_dS2=d_B_dS2*betahat;
    dYnew_dS2=dS2 - dfB_dS2;
    dfB_dS2pde=d_B_dS2*beta_pde;
    dYnew_dS2pde=dS2 - dfB_dS2pde;
    
    %% ESK surface
    Ynew=Y - B1*betahat;
    f_ESK(:,i)=Wm_price.ESK'*[Ynew;dYnew_dS1;dYnew_dS2]+Btest*Wm_beta.ESK;
    %% ESKpde surface
    Ynew_pde=Y - B1*beta_pde;
    f_ESK_pde(:,i)=Wm_price.ESKpde'*[Ynew_pde;dYnew_dS1pde;dYnew_dS2pde]+Btest*Wm_beta.ESKpde

%% Greeks
%% ESK
    delta1_ESK(:,i) = Wm_delta1.ESK'*[Ynew;dYnew_dS1;dYnew_dS2] + d_Btest_delta*betahat;
    delta2_ESK(:,i) = Wm_delta2.ESK'*[Ynew;dYnew_dS1;dYnew_dS2] + d_Btest_delta*betahat;
    gamma11_ESK(:,i) = Wm_gamma11.ESK'*[Ynew;dYnew_dS1;dYnew_dS2] + d_Btest_gamma*betahat;
    gamma12_ESK(:,i) = Wm_gamma12.ESK'*[Ynew;dYnew_dS1;dYnew_dS2] + d_Btest_gamma*betahat;
    gamma21_ESK(:,i) = Wm_gamma12.ESK'*[Ynew;dYnew_dS1;dYnew_dS2] + d_Btest_gamma*betahat;
    gamma22_ESK(:,i) = Wm_gamma22.ESK'*[Ynew;dYnew_dS1;dYnew_dS2] + d_Btest_gamma*betahat;
    theta_ESK(:,i) = Wm_theta.ESK'*[Ynew;dYnew_dS1;dYnew_dS2] + d_Btest_theta*betahat;

%% ESKpde
    delta1_ESK_pde(:,i) = Wm_delta1.ESKpde'*[Ynew_pde;dYnew_dS1pde;dYnew_dS2pde] + d_Btest_delta*Wm_beta.ESKpde;
    delta2_ESK_pde(:,i) = Wm_delta2.ESKpde'*[Ynew_pde;dYnew_dS1pde;dYnew_dS2pde] + d_Btest_delta*Wm_beta.ESKpde;
    gamma11_ESK_pde(:,i) = Wm_gamma11.ESKpde'*[Ynew_pde;dYnew_dS1pde;dYnew_dS2pde] + d_Btest_gamma*Wm_beta.ESKpde;
    gamma12_ESK_pde(:,i) = Wm_gamma12.ESKpde'*[Ynew_pde;dYnew_dS1pde;dYnew_dS2pde] + d_Btest_gamma*Wm_beta.ESKpde;
    gamma21_ESK_pde(:,i) = Wm_gamma12.ESKpde'*[Ynew_pde;dYnew_dS1pde;dYnew_dS2pde] + d_Btest_gamma*Wm_beta.ESKpde;
    gamma22_ESK_pde(:,i) = Wm_gamma22.ESKpde'*[Ynew_pde;dYnew_dS1pde;dYnew_dS2pde] + d_Btest_gamma*Wm_beta.ESKpde;
    theta_ESK_pde(:,i) = Wm_theta.ESKpde'*[Ynew_pde;dYnew_dS1pde;dYnew_dS2pde] + d_Btest_theta*Wm_beta.ESKpde;
end

Port_true=sum(Y_true,2);Port_T_true=sum(theta_true,2);
Port_ESK=sum(f_ESK,2);Port_ESK_pde=sum(f_ESK_pde,2);
Port_T_ESK=sum(theta_ESK,2);Port_T_ESK_pde=sum(theta_ESK_pde,2);


RMSE_V.ESK = sqrt(mean((Port_ESK - Port_true).^2));
RMSE_V.ESKpde = sqrt(mean((Port_ESK_pde - Port_true).^2));


RMSE_theta.ESK = sqrt(mean((Port_T_ESK - Port_T_true).^2));
RMSE_theta.ESKpde = sqrt(mean((Port_T_ESK_pde - Port_T_true).^2));

RMSE_delta1.ESK=zeros(40,1);RMSE_delta2.ESK=zeros(40,1);
RMSE_delta1.ESKpde=zeros(40,1);RMSE_delta2.ESKpde=zeros(40,1);
RMSE_gamma11.ESK=zeros(40,1);RMSE_gamma11.ESKpde=zeros(40,1);
RMSE_gamma12.ESK=zeros(40,1);RMSE_gamma12.ESKpde=zeros(40,1);
RMSE_gamma21.ESK=zeros(40,1);RMSE_gamma21.ESKpde=zeros(40,1);
RMSE_gamma22.ESK=zeros(40,1);RMSE_gamma22.ESKpde=zeros(40,1);
for i=1:40
    RMSE_delta1.ESK(i)=sqrt(mean((delta1_ESK(:,i) - delta_true(:,1,i)).^2));
    RMSE_delta2.ESK(i)=sqrt(mean((delta2_ESK(:,i) - delta_true(:,2,i)).^2));
    RMSE_delta1.ESKpde(i)=sqrt(mean((delta1_ESK_pde(:,i) - delta_true(:,1,i)).^2));
    RMSE_delta2.ESKpde(i)=sqrt(mean((delta2_ESK_pde(:,i) - delta_true(:,2,i)).^2));
    RMSE_gamma11.ESK(i)=sqrt(mean((gamma11_ESK(:,i) - gamma_true(:,1,i)).^2));
    RMSE_gamma12.ESK(i)=sqrt(mean((gamma12_ESK(:,i) - gamma_true(:,2,i)).^2));
    RMSE_gamma21.ESK(i)=sqrt(mean((gamma21_ESK(:,i) - gamma_true(:,3,i)).^2));
    RMSE_gamma22.ESK(i)=sqrt(mean((gamma22_ESK(:,i) - gamma_true(:,4,i)).^2));
    RMSE_gamma11.ESKpde(i)=sqrt(mean((gamma11_ESK_pde(:,i) - gamma_true(:,1,i)).^2));
    RMSE_gamma12.ESKpde(i)=sqrt(mean((gamma12_ESK_pde(:,i) - gamma_true(:,2,i)).^2));
    RMSE_gamma21.ESKpde(i)=sqrt(mean((gamma21_ESK_pde(:,i) - gamma_true(:,3,i)).^2));
    RMSE_gamma22.ESKpde(i)=sqrt(mean((gamma22_ESK_pde(:,i) - gamma_true(:,4,i)).^2));  
end

ARMSE_delta1.ESK=mean(RMSE_delta1.ESK);ARMSE_delta2.ESK=mean(RMSE_delta2.ESK);
ARMSE_delta1.ESKpde=mean(RMSE_delta1.ESKpde);ARMSE_delta2.ESKpde=mean(RMSE_delta2.ESKpde);
ARMSE_gamma11.ESK=mean(RMSE_gamma11.ESK);ARMSE_gamma11.ESKpde=mean(RMSE_gamma11.ESKpde);
ARMSE_gamma12.ESK=mean(RMSE_gamma12.ESK);ARMSE_gamma12.ESKpde=mean(RMSE_gamma12.ESKpde);
ARMSE_gamma21.ESK=mean(RMSE_gamma21.ESK);ARMSE_gamma21.ESKpde=mean(RMSE_gamma21.ESKpde);
ARMSE_gamma22.ESK=mean(RMSE_gamma22.ESK);ARMSE_gamma22.ESKpde=mean(RMSE_gamma22.ESKpde);





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
    
    

