%% consider the portfolio has one Basket option
function [Wm_price,Wm_delta1,Wm_delta2,Wm_gamma11,Wm_gamma12,Wm_gamma22,Wm_theta,Wm_beta,finpar]=basket_kriging_PDE_grad(r,sigma,S0,K,T,Rho,Nd,X,Y,vY,Vmatrix,Xtest,seed,poly_d,dS1,dS2)
%tic
%% generate testing points 
p=sobolset(3);
%% generate approximating points for pdf-objective function
S_taumind=S0*0.48;S_taumaxd=S0*1.52;
t_mind=0;t_maxd=1;
u2=p(502:501+Nd,:);
Sd = S_taumind+(S_taumaxd-S_taumind).*u2(:,1:2);
td=t_mind+(t_maxd-t_mind)*u2(:,3);
Xd=[Sd,td];
B1=ones(size(X,1),1);
gammaP=2;
para_pen=0.001;%penality parameter in solving MLE
SKmodel1=SKfit(X, Y, B1, vY, gammaP,para_pen);% original SK

SKpar=[SKmodel1.beta;SKmodel1.tausquared;SKmodel1.theta];

SKmodel_pde=SKfit_pde_basket2(X, Y, B1, vY, gammaP,para_pen,Xd,r,sigma,Rho,dS1,dS2,Vmatrix);%pde enhanced SK

SKpdepar=[SKmodel_pde.beta;SKmodel_pde.tausquared;SKmodel_pde.theta];

finpar=[SKpar,SKpdepar];

Btest=polybasis(Xtest,poly_d);

[w,betahat] = SKpredict_M_linear(SKmodel1,Xtest,Btest);%betahat=0;
[w_pde,betahat_pde] = SKpredict_M_linear(SKmodel_pde,Xtest,Btest);
[wc_pde,wc_delta1_pde,wc_delta2_pde,wc_gamma11_pde,wc_gamma12_pde,wc_gamma22_pde,wc_theta_pde,betahatc_pde] = SKpredict_cokriging_M_linear_delta_pde2(SKmodel_pde,Xtest,Btest,Vmatrix);
[wc_e,wc_delta1_e,wc_delta2_e,wc_gamma11_e,wc_gamma12_e,wc_gamma22_e,wc_theta_e,betahatc_e] = SKpredict_cokriging_M_linear_delta_pde2(SKmodel1,Xtest,Btest,Vmatrix);

Wm_price.ESKpde=wc_pde;Wm_price.ESK=wc_e;
Wm_theta.ESKpde=wc_theta_pde;Wm_theta.ESK=wc_theta_e;
Wm_delta1.ESKpde=wc_delta1_pde;Wm_delta1.ESK=wc_delta1_e;
Wm_delta2.ESKpde=wc_delta2_pde;Wm_delta2.ESK=wc_delta2_e;
Wm_gamma11.ESKpde=wc_gamma11_pde;Wm_gamma11.ESK=wc_gamma11_e;
Wm_gamma22.ESKpde=wc_gamma22_pde;Wm_gamma22.ESK=wc_gamma22_e;
Wm_gamma12.ESKpde=wc_gamma12_pde;Wm_gamma12.ESK=wc_gamma12_e;
Wm_beta.ESKpde=betahatc_pde;Wm_beta.ESK=betahatc_e;
Wm_beta.SK=betahat;Wm_beta.SKpde=betahat_pde;



    
    

