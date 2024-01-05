%% consider the portfolio has one European option
%% in this example, we test the mutually cokriging methods
function [V,VdS0]=KrigEuro_cokriging_delta_Asian1(X,K,T,iday,n,r_n,sigma_vg,nu_vg,theta_vg,M,poly_d,Xtest)
% alpha is the range precentage of the S0
% n is # of design points
addpath ./SK_revised2

%[Y,vY,dY,vdY,vYdY] = VG_call_Asian_v_hedge(X(:,1),K,T,iday,n,M,r_n,X(:,2),nu_vg,X(:,3));
[Y,vY,dY,vdY,vYdY] = VG_call_Asian_v_hedge(X(:,1),K,T,iday,n,M,r_n,sigma_vg,nu_vg,theta_vg);

vMatrix = [diag(vY),diag(vYdY);diag(vYdY),diag(vdY)];

B1=polybasis(X,poly_d);
gammaP=2;
SKmodel1=SKfit(X, Y, B1, vY, gammaP,0.5);
SKmodel1_der=SKfit_der(X, dY, B1, vdY, gammaP,0.1);
SKmodel_enhanced=SK_grad1_fit(X, Y, dY, B1, vY, vMatrix, gammaP,0.1);

%Xtest=linspace(min(X),max(X),100)';


Btest=polybasis(Xtest,poly_d);
d_B=d_polybasis(X,poly_d);
d_B_only=polybasis(X,poly_d);
d_Btest=d_polybasis(Xtest,poly_d);
d_Btest_only=polybasis(Xtest,poly_d);
[w,betahat] = SKpredict_M_linear(SKmodel1,Xtest,Btest);
[w_der,betahat_der] = SKpredict_M_linear_der(SKmodel1_der,Xtest,d_Btest_only);

[w_enhanced,w_d_enhanced,~] = SKpredict_cokriging_M_linear(SKmodel_enhanced,Xtest,Btest);


dfB=d_B*betahat;
dYnew=dY - dfB;
Ynew=Y - B1*betahat;
f_SK=w'*Ynew+Btest*betahat;
f_ESK=w_enhanced'*[Ynew;dYnew]+Btest*betahat;
df_ESK=w_d_enhanced'*[Ynew;dYnew]+d_Btest*betahat;

dY_onlynew = dY - d_B_only*betahat_der;
df_SK = w_der'*dY_onlynew + d_Btest_only*betahat_der;

V.SK = f_SK; V.ESK = f_ESK; VdS0.SK = df_SK; VdS0.ESK = df_ESK;
% 
% subplot(1,2,1)
% plot(X,Y,'k*',Xtest,f_SK,'b--',Xtest,f_ESK,'r-');
% subplot(1,2,2)
% plot(X,dY,'k*',Xtest,df_SK,'b--',Xtest,df_ESK,'r-');








%%hedging cases
% X0=Xtest(2:end);
% dX0=Xtest(1:end-1)-X0;
% dphi_true = (Y_true(2:end) - Y_true(1:end-1));
% dphi_enhanced = (f_enhanced(2:end) - f_enhanced(1:end-1));
% dphi_only = (f2(2:end) - f2(1:end-1));
% dphi_true_hedge = dY_true(1:end-1).*dX0;
% dphi_enhanced_hedge = df_enhanced(1:end-1).*dX0;
% dphi_only_hedge = df_only(1:end-1).*dX0;
% figure(2)
% subplot(1,2,1)
% plot(X0,dphi_true-dphi_true_hedge,'k-',X0,dphi_true-dphi_enhanced_hedge,'g--',X0,dphi_true-dphi_only_hedge,'b--','LineWidth',2)
% legend('Real Gap of real \Delta hedged','Real Gap of enhanced \Delta hedged','Real Gap of only \Delta hedged')
% subplot(1,2,2)
% plot(X0,dphi_true-dphi_true_hedge,'k-',X0,dphi_enhanced-dphi_enhanced_hedge,'g--',X0,dphi_only-dphi_only_hedge,'b--','LineWidth',2)
% legend('Real Gap of real \Delta hedged','Model Gap of enhanced \Delta hedged','Model Gap of only \Delta hedged')




function basis=polybasis(X,deg)
basis=ones(size(X,1),1);
for i=1:deg
    basis=[basis,X.^i];
end

function d_basis=d_polybasis(X,deg)
if deg==0
    d_basis=0*X;
else
d_basis=[];
for i=1:deg+1
    d_basis=[d_basis,(i-1)*X.^(i-2)];
end
end


    
    

