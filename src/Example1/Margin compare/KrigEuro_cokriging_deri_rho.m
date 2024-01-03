%% consider the portfolio has one European option
%% in this example, we test the mutually cokriging methods
function []=KrigEuro_cokriging_deri_rho(S0,r,sigma,T,K,n,alpha,M,t,poly_d,Tru_MC_idx)
% alpha is the range precentage of the S0
% n is # of design points
addpath ./SK_revised2

rng(1)
%X=linspace(1-alpha,1+alpha,n)'*S0;
%Xtem1=linspace(-alpha,alpha,n);sign1=ones(1,n);sign1(1:ceil(n/2))=-1;X=(1+1/alpha*Xtem1.^2.*sign1)'*S0;
%X=[97;101;103];
%X=[0.001;0.011;0.022;0.031;0.039;0.049];
%X=[0.004;0.021;0.042;0.061;0.080;0.095];
u1=lhsdesign(n,1);
X=0.001+u1*0.099;
Xtest=(0.001:0.001:0.1)';
%Xtest=(95:0.5:105)';
%Xtest=(98:0.5:102)';
n=length(X);
if Tru_MC_idx==1
    [Y,dY]=bls_Euro_call_rho(S0,K,X,sigma,t,T);
    vY=0.01*ones(n,1);
    vdY=0.01*ones(n,1);
    vMatrix = 0.01*diag(ones(2*n,1));%%%make sure the inverse exists
else
    [Y,vY,dY,vdY,vYdY]=EuroMC_rho(S0,X,sigma,T,K,M);
    vMatrix = [diag(vY),diag(vYdY);diag(vYdY),diag(vdY)];
end

%[Yt,dYt]=bls_Euro_call_rho(S0,K,X,sigma,t,T);
B1=polybasis(X,poly_d);
gammaP=2;
SKmodel1=SKfit(X, Y, B1, vY, gammaP,0.1);
SKmodel1_der=SKfit_der(X, dY, B1, vdY, gammaP,0.1);
SKmodel_enhanced=SK_grad1_fit(X, Y, dY, B1, vY, vMatrix, gammaP,0.1);

[Y_true,dY_true]=bls_Euro_call_rho(S0,K,Xtest,sigma,t,T);
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
f2=w'*Ynew+Btest*betahat;
f_enhanced=w_enhanced'*[Ynew;dYnew]+Btest*betahat;
%df_co=w_d_enhanced'*[Ynew;dYnew]+d_Btest*betahat;


dY_onlynew = dY - d_B_only*betahat_der;
df_only = w_der'*dY_onlynew + d_Btest_only*betahat_der;

idx=1;
dw = SKpredict_M_linear_dw(SKmodel1,Xtest,Btest,idx);

%dw_enhanced = SKpredict_grad1_M_linear_dw(SKmodel_enhanced,Xtest,Btest,idx);

df2=dw'*Ynew + d_Btest*betahat;

df_enhanced = w_d_enhanced' * [Ynew;dYnew] + d_Btest*betahat;

figure(1)
%plot(X,Y,'kx',Xtest,Y_true,'r',Xtest,f1,'b--',Xtest,f2,'c:','LineWidth',2)
%legend('point','True','Barry','Linearcomb')
plot(X,Y,'kx',Xtest,Y_true,'k',Xtest,f_enhanced,'r--',Xtest,f2,'b:','LineWidth',2,'MarkerSize',12)
legend('point','True','GESK','SK')
set(gca,'FontSize',16);
xlabel('r')
ylabel('Price')
figure(2)
plot(X,dY,'kx',Xtest,dY_true,'k',Xtest,df_enhanced,'r--',Xtest,df_only,'b:',Xtest, df2,'c:','LineWidth',2,'MarkerSize',12)
legend('point','True','GESK','SK','SK-D')
set(gca,'FontSize',16);
xlabel('r')
ylabel('Rho')


%%hedging cases
% X0=Xtest(1:end-1);
% dX0=Xtest(2:end)-X0;
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


function [V_c,dV_c]=bls_Euro_call_rho(S,K,r,sigma,t,T)
tau=T-t;
d1=1./(sigma*sqrt(tau)).*(log(S/K)+(r+sigma.^2/2)*tau);
d2=d1-sigma*sqrt(tau);
V_c = normcdf(d1,0,1).*S-normcdf(d2,0,1).*K.*exp(-r*tau);
dV_c= K*tau*exp(-r*tau).*normcdf(d2,0,1);





function [V,vV,G,vG,vVG]=EuroMC_rho(S0,r,sigma,T,K,M)
M = ones(size(r))*M;%% larger sigma has larger variance
%M=round(sigma*10*M);
l=length(r);
V=zeros(l,1);vV=zeros(l,1);G=zeros(l,1);vG=zeros(l,1);vVG=zeros(l,1);
for i=1:l
Z=normrnd(0,1,M(i),1);
ST=S0*exp((r(i)-0.5*sigma^2)*T+sigma*sqrt(T)*Z);
Payoff=exp(-r(i)*T)*max(ST-K,0);
V(i)=mean(Payoff);
vV(i)=var(Payoff)/M(i);
IPA_rho=-T*exp(-r(i)*T)*(ST - K).*(ST>K)+ exp(-r(i)*T)*T*ST.*(ST>K);
G(i)=mean(IPA_rho);
vG(i)=var(IPA_rho)/M(i);
covVG = cov(IPA_rho,Payoff);
vVG(i) = covVG(1,2)/M(i);
end



function basis=polybasis(X,deg)
basis=[];
for i=1:deg+1
    basis=[basis,X.^(i-1)];
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

    
    

