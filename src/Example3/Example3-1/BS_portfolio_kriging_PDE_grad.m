%% consider the portfolio has one Basket option
function [Wm_price,Wm_delta,Wm_gamma,Wm_theta,Wm_beta,finpar]=BS_portfolio_kriging_PDE_grad(r,sigma,S0,T,Nd,X,Y,vY,Vmatrix,Xtest,poly_d,dS)
%tic
%% generate testing points 
p=sobolset(2);
%u=p(2:Nt+1,:);
%Stest = ones(Nt,1)*S_taumin+(ones(Nt,1)*S_taumax-ones(Nt,1)*S_taumin).*u(:,1:2);
%ttest=t_min+(t_max-t_min).*u(:,3);
%Xtest=[Stest,ttest];
%[Y_true,delta_true,gamma_true,theta_true]=portfolio_basket_price(Xtest(:,1:2),r,sigma,T-Xtest(:,3),rho,Nt,K1,K2,K3);
%Y_true is Nt*1, delta_true is Nt*2, gamma_true is Nt*4 (d1d1,d1d2,d2d1,d2d2), theta_true is Nt*1
%[Y1_true,delta1_true,gamma1_true,theta1_true]=portfolio_basket_price2(Xtest(:,1:2),r,sigma,T-Xtest(:,3),Rho,K1,K2,K3);
%[Y2_true,delta2_true,gamma2_true,theta2_true]=portfolio_basket_price2(Xtest(:,1:2),r,sigma,T-Xtest(:,3),Rho,K4,K5,K6);
%Y_true=10*Y1_true+Y2_true;delta_true=10*delta1_true+delta2_true;gamma_true=10*gamma1_true+gamma2_true;theta_true=-10*theta1_true-theta2_true;
%% generate approximating points for pdf-objective function
S_taumind=S0*0.8;S_taumaxd=S0*1.2;
t_mind=0;t_maxd=T;
u2=p(502:501+Nd,:);
Sd = S_taumind+(S_taumaxd-S_taumind).*u2(:,1);
td=t_mind+(t_maxd-t_mind)*u2(:,2);
Xd=[Sd,td];

%% Generate design points

%B1=polybasis(X,poly_d);
B1=ones(size(X,1),1);
gammaP=2;
para_pen=0.0;%penality parameter in solving MLE
SKmodel1=SKfit(X, Y, B1, vY, gammaP,para_pen);% original SK

SKpar=[SKmodel1.beta;SKmodel1.tausquared;SKmodel1.theta];

%xlswrite('Para.xls',SKmodel1.beta,seed,'A1')
%xlswrite('Para.xls',SKmodel1.tausquared,seed,'A2')
%xlswrite('Para.xls',SKmodel1.theta,seed,'A3:A5')

SKmodel_pde=SKfit_pde(X, Y, B1, vY, gammaP,para_pen,Xd,r,sigma,dS,Vmatrix);%pde enhanced SK

SKpdepar=[SKmodel_pde.beta;SKmodel_pde.tausquared;SKmodel_pde.theta];

finpar=[SKpar,SKpdepar];
Btest=polybasis(Xtest,poly_d);
d_Btest_only=polybasis(Xtest,poly_d);
[w,betahat] = SKpredict_M_linear(SKmodel1,Xtest,Btest);%betahat=0;
[w_pde,betahat_pde] = SKpredict_M_linear(SKmodel_pde,Xtest,Btest);%beta_pde=0;
[w_delta,w_gamma,w_theta,~] = SKpredict_M_linear_der_delta_gamma_theta_pde(SKmodel1,Xtest,d_Btest_only);
[w_delta_pde,w_gamma_pde,w_theta_pde,~]=SKpredict_M_linear_der_delta_gamma_theta_pde(SKmodel_pde,Xtest,Btest);
[wc_pde,wc_delta_pde,wc_gamma_pde,wc_theta_pde,betahatc_pde] = SKpredict_cokriging_M_linear_delta_pde2(SKmodel_pde,Xtest,Btest,Vmatrix);
[wc_e,wc_delta_e,wc_gamma_e,wc_theta_e,betahatc_e] = SKpredict_cokriging_M_linear_delta_pde2(SKmodel1,Xtest,Btest,Vmatrix);

Wm_price.SK=w;Wm_price.SKpde=w_pde;Wm_price.ESKpde=wc_pde;Wm_price.ESK=wc_e;
Wm_theta.SK=w_theta;Wm_theta.SKpde=w_theta_pde;Wm_theta.ESKpde=wc_theta_pde;Wm_theta.ESK=wc_theta_e;
Wm_delta.SK=w_delta;Wm_delta.SKpde=w_delta_pde;Wm_delta.ESKpde=wc_delta_pde;Wm_delta.ESK=wc_delta_e;
Wm_gamma.SK=w_gamma;Wm_gamma.SKpde=w_gamma_pde;Wm_gamma.ESKpde=wc_gamma_pde;Wm_gamma.ESK=wc_gamma_e;
Wm_beta.SK=betahat;Wm_beta.SKpde=betahat_pde;Wm_beta.ESKpde=betahatc_pde;Wm_beta.ESK=betahatc_e;



    
    

