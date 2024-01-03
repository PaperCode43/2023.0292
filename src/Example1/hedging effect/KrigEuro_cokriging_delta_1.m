%% consider the portfolio has one European option
%% in this example, we test the mutually cokriging methods
function [Y_true,dY_true,f_enhanced,df_enhanced,f_SK,df_SK]=KrigEuro_cokriging_delta_1(X,r,sigma,T,K,M,t,poly_d,Tru_MC_idx,Xtest)
% alpha is the range precentage of the S0
% n is # of design points
addpath ./SK_revised2

n=length(X);
if Tru_MC_idx==1
    [Y,dY]=bls_Euro_call_delta(X,K,r,sigma,t,T);
    vY=0.0*ones(n,1);
    vdY=0.000*ones(n,1);
    vMatrix = diag([0*ones(n,1);0.000*ones(n,1)]);%%%make sure the inverse exists
else
    [Y,vY,dY,vdY,vYdY]=EuroMC(X,r,sigma,T,K,M);
    vMatrix = [diag(vY),diag(vYdY);diag(vYdY),diag(vdY)]+0*diag(ones(2*n,1));
end


B1=polybasis(X,poly_d);
gammaP=2;
SKmodel1=SKfit(X, Y, B1, vY, gammaP,0.1);
SKmodel1_der=SKfit_der(X, dY, B1, vdY, gammaP,0.1);
SKmodel_enhanced=SK_grad1_fit(X, Y, dY, B1, vY, vMatrix, gammaP,0.1);

[Y_true,dY_true]=bls_Euro_call_delta(Xtest,K,r,sigma,t,T);
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
f_enhanced=w_enhanced'*[Ynew;dYnew]+Btest*betahat;
df_enhanced=w_d_enhanced'*[Ynew;dYnew]+d_Btest*betahat;


dY_onlynew = dY - d_B_only*betahat_der;
df_SK = w_der'*dY_onlynew + d_Btest_only*betahat_der;



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


function [V_c,dV_c]=bls_Euro_call_delta(S,K,r,sigma,t,T)
tau=T-t;
d1=1/(sigma*sqrt(tau))*(log(S./K)+(r+sigma^2/2)*tau);
d2=d1-sigma*sqrt(tau);
V_c = normcdf(d1,0,1).*S-normcdf(d2,0,1).*K*exp(-r*tau);
dV_c= normcdf(d1,0,1);




function [V,vV,G,vG,vVG]=EuroMC(S0,r,sigma,T,K,M)
l=length(S0);
V=zeros(l,1);vV=zeros(l,1);G=zeros(l,1);vG=zeros(l,1);vVG=zeros(l,1);
for i=1:l
ST=S0(i)*exp((r-0.5*sigma^2)*T+sigma*sqrt(T)*normrnd(0,1,M,1));
Payoff=exp(-r*T)*max(ST-K,0);
V(i)=mean(Payoff);
vV(i)=var(Payoff)/M;
IPA_S0=exp(-r*T)*ST/S0(i).*(ST>K);
G(i)=mean(IPA_S0);
vG(i)=var(IPA_S0)/M;
covVG = cov(IPA_S0,Payoff);
vVG(i) = covVG(1,2)/M;
end



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

    
    

