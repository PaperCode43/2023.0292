%% consider the portfolio has one European option
%% in this example, we test the mutually cokriging methods
function [RMSE_V,RMSE_dS0,RMSE_dsigma,RMSE_dtheta]=KrigEuro_cokriging_Lookback(X,K,T,r_n,nu_vg,M,poly_d,Xtest,Ytest_true)
% alpha is the range precentage of the S0
% n is # of design points
addpath ./SK_revised2
addpath ./datafiles

n=size(X,1);

[Y,vY,dS0,vdS0,dsigma,vdsigma,dtheta,vdtheta,vMatrix] = VG_call_Lookback_v(X(:,1),K,T,n,M,r_n,X(:,2),nu_vg,X(:,3));


B1=polybasis(X,poly_d);
gammaP=2;
SKmodel1=SKfit(X, Y, B1, vY, gammaP,0.1);
SKmodel1_dS0=SKfit_der(X, dS0, B1, vdS0, gammaP,0.1);
SKmodel1_dsigma=SKfit_der(X, dsigma, B1, vdsigma, gammaP,0.1);
SKmodel1_dtheta=SKfit_der(X, dtheta, B1, vdtheta, gammaP,0.1);
SKmodel_enhanced=SK_grad1_fit_allgreek_Asian(X, Y, dS0, B1, vY, vMatrix, gammaP,0.1);%%delta is a invalid input see code of SK_grad1_fit
                                                                                 %%Also, this command is to learn parameters, also just use Y
Btest=polybasis(Xtest,poly_d);
d_Btest_only=polybasis(Xtest,poly_d);
[w,betahat] = SKpredict_M_linear(SKmodel1,Xtest,Btest);

[w_dS0,betahat_dS0] = SKpredict_M_linear_der(SKmodel1_dS0,Xtest,d_Btest_only);
[w_dsigma,betahat_dsigma] = SKpredict_M_linear_der(SKmodel1_dsigma,Xtest,d_Btest_only);
[w_dtheta,betahat_dtheta] = SKpredict_M_linear_der(SKmodel1_dtheta,Xtest,d_Btest_only);

[w_enhanced,w_delta_enhanced,w_vega_enhanced,w_rho_enhanced,~]...
    = SKpredict_cokriging_M_linear_allgreek_Asian(SKmodel_enhanced,Xtest,Btest);

%%% SK surface 
Ynew=Y - B1*betahat;
f_SK=w'*Ynew+Btest*betahat;
%%% dS0
d_B_dS0=d_polybasis(X,poly_d,1);
d_Btest_dS0=d_polybasis(Xtest,poly_d,1);
dfB_dS0=d_B_dS0*betahat_dS0;
dYnew_dS0=dS0 - dfB_dS0;
dS0_SK = w_dS0'*dYnew_dS0 + d_Btest_dS0*betahat_dS0;
%%% dsigma
d_B_dsigma=d_polybasis(X,poly_d,2);
d_Btest_dsigma=d_polybasis(Xtest,poly_d,2);
dfB_dsigma=d_B_dsigma*betahat_dsigma;
dYnew_dsigma=dsigma - dfB_dsigma;
dsigma_SK = w_dsigma'*dYnew_dsigma + d_Btest_dsigma*betahat_dsigma;
%%% rho
d_B_dtheta=d_polybasis(X,poly_d,3);
d_Btest_dtheta=d_polybasis(Xtest,poly_d,3);
dfB_dtheta=d_B_dtheta*betahat_dtheta;
dYnew_dtheta=dtheta - dfB_dtheta;
dtheta_SK = w_dtheta'*dYnew_dtheta + d_Btest_dtheta*betahat_dtheta;

f_enhanced=w_enhanced'*[Ynew;dYnew_dS0;dYnew_dsigma;dYnew_dtheta]+Btest*betahat;
df_dS0_enhanced=w_delta_enhanced'*[Ynew;dYnew_dS0;dYnew_dsigma;dYnew_dtheta]+d_Btest_dS0*betahat;
df_dsigma_enhanced=w_vega_enhanced'*[Ynew;dYnew_dS0;dYnew_dsigma;dYnew_dtheta]+d_Btest_dsigma*betahat;
df_dtheta_enhanced=w_rho_enhanced'*[Ynew;dYnew_dS0;dYnew_dsigma;dYnew_dtheta]+d_Btest_dtheta*betahat;

Y_true=Ytest_true.V;
dS0_true = Ytest_true.dS0;
dsigma_true = Ytest_true.dsigma;
dtheta_true = Ytest_true.dtheta;

RMSE_V.SK = sqrt(mean((f_SK - Y_true).^2));
RMSE_V.ESK = sqrt(mean((f_enhanced - Y_true).^2));

RMSE_dS0.SK = sqrt(mean((dS0_SK - dS0_true).^2));
RMSE_dS0.ESK = sqrt(mean((df_dS0_enhanced - dS0_true).^2));

RMSE_dsigma.SK = sqrt(mean((dsigma_SK - dsigma_true).^2));
RMSE_dsigma.ESK = sqrt(mean((df_dsigma_enhanced - dsigma_true).^2));

RMSE_dtheta.SK = sqrt(mean((dtheta_SK - dtheta_true).^2));
RMSE_dtheta.ESK = sqrt(mean((df_dtheta_enhanced - dtheta_true).^2));




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
    
    

