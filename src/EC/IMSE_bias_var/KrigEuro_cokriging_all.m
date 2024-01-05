%% consider the portfolio has one European option
%% in this example, we test the mutually cokriging methods
function [RMSE_V,Bias_V,std_V,RMSE_delta,Bias_delta,std_delta,RMSE_vega,Bias_vega,std_vega,RMSE_rho,Bias_rho,std_rho,RMSE_theta,Bias_theta,std_theta]=KrigEuro_cokriging_all(X,K,M,t,poly_d,Tru_MC_idx)
% alpha is the range precentage of the S0
% n is # of design points
addpath ./SK_revised2
%%%% S0 [80 120], sigma [0.01 0.3], r [0.001,0.1], T [0.01,2] 

%% generate testing points 
p=sobolset(4);
u=p(2:1001,:);
%%%%X=(S0,sigma,r,T)
Xtest=zeros(1000,4);
Xtest(:,1) = (1+(2*u(:,1)-1)*0.2)*100;
Xtest(:,2) = u(:,2)*0.3;
Xtest(:,3) = u(:,3)*0.1;
Xtest(:,4) = u(:,4)*2;
%Xtest=X;



%Xtest=[99.2,98.8,97.6,97.4,98.1,98.9,99.4,100.1,101.0,101.9,102.4,103.1,103.9,104.5]';
%Xtest=(98:0.5:102)';
n=size(X,1);
Y=zeros(n,1);delta=zeros(n,1);rho=zeros(n,1);theta=zeros(n,1);vega=zeros(n,1);
if Tru_MC_idx==1
    for i=1:n
    [Y(i),delta(i),vega(i),rho(i),theta(i)]=bls_Euro_call_greek(X(i,1),K,X(i,3),X(i,2),t,X(i,4));
    end   
    vY=0.0*ones(n,1);
    vdelta=0.000*ones(n,1);
    vvega=0.000*ones(n,1);
    vrho=0.000*ones(n,1);
    vtheta=0.000*ones(n,1);
    vMatrix = diag([0*ones(n,1);0*ones(n,1);0*ones(n,1);0*ones(n,1);0*ones(n,1)]);%%%make sure the inverse exists
else
    [Y,vY,delta,vdelta,vega,vvega,rho,vrho,theta,vtheta,vMatrix]=EuroMC_greek(X(:,1),X(:,3),X(:,2),X(:,4),K,M);
    %[Y,vY,dY,vdY,vYdY]=EuroMC(X,r,sigma,T,K,M);
    %vMatrix = [diag(vY),diag(vYdY);diag(vYdY),diag(vdY)]+0*diag(ones(2*n,1));
end


B1=polybasis(X,poly_d);
gammaP=2;
fprintf('SK begin \n');
SKmodel1=SKfit(X, Y, B1, vY, gammaP,0.1);
fprintf('SK delta begin \n');
SKmodel1_delta=SKfit_der(X, delta, B1, vdelta, gammaP,0.1);
fprintf('SK vega begin \n');
SKmodel1_vega=SKfit_der(X, vega, B1, vvega, gammaP,0.1);
fprintf('SK rho begin \n');
SKmodel1_rho=SKfit_der(X, rho, B1, vrho, gammaP,0.1);
fprintf('SK theta begin \n');
SKmodel1_theta=SKfit_der(X, theta, B1, vtheta, gammaP,0.1);
fprintf('SK all begin \n');
SKmodel_enhanced=SK_grad1_fit_allgreek(X, Y, delta, B1, vY, vMatrix, gammaP,0.1);%%delta is a invalid input see code of SK_grad1_fit
fprintf('SK allconstant begin \n');
SKmodel_enhanced_const=SK_grad1_fit_allgreek_const(X, Y, delta, B1, vY, vMatrix, gammaP,0.1);

[Y_true,delta_true,vega_true,rho_true,theta_true]=bls_Euro_call_greek(Xtest(:,1),K,Xtest(:,3),Xtest(:,2),t,Xtest(:,4));

Btest=polybasis(Xtest,poly_d);
%d_B=d_polybasis(X,poly_d);
%d_B_only=polybasis(X,poly_d);
%d_Btest=d_polybasis(Xtest,poly_d);
d_Btest_only=polybasis(Xtest,poly_d);
[w,betahat] = SKpredict_M_linear(SKmodel1,Xtest,Btest);

[w_delta,betahat_delta] = SKpredict_M_linear_der(SKmodel1_delta,Xtest,d_Btest_only);
[w_vega,betahat_vega] = SKpredict_M_linear_der(SKmodel1_vega,Xtest,d_Btest_only);
[w_rho,betahat_rho] = SKpredict_M_linear_der(SKmodel1_rho,Xtest,d_Btest_only);
[w_theta,betahat_theta] = SKpredict_M_linear_der(SKmodel1_theta,Xtest,d_Btest_only);

%[w_enhanced,w_d_enhanced,~] = SKpredict_cokriging_M_linear(SKmodel_enhanced,Xtest,Btest);

[w_enhanced,w_delta_enhanced,w_vega_enhanced,w_rho_enhanced,w_theta_enhanced,~]...
    = SKpredict_cokriging_M_linear_allgreek(SKmodel_enhanced,Xtest,Btest);

%%% SK surface 
Ynew=Y - B1*betahat;
f_SK=w'*Ynew+Btest*betahat;
%%% delta
d_B_delta=d_polybasis(X,poly_d,1);
d_Btest_delta=d_polybasis(Xtest,poly_d,1);
dfB_delta=d_B_delta*betahat_delta;
dYnew_delta=delta - dfB_delta;
delta_SK = w_delta'*dYnew_delta + d_Btest_delta*betahat_delta;
%%% vega
d_B_vega=d_polybasis(X,poly_d,2);
d_Btest_vega=d_polybasis(Xtest,poly_d,2);
dfB_vega=d_B_vega*betahat_vega;
dYnew_vega=vega - dfB_vega;
vega_SK = w_vega'*dYnew_vega + d_Btest_vega*betahat_vega;
%%% rho
d_B_rho=d_polybasis(X,poly_d,3);
d_Btest_rho=d_polybasis(Xtest,poly_d,3);
dfB_rho=d_B_rho*betahat_rho;
dYnew_rho=rho - dfB_rho;
rho_SK = w_rho'*dYnew_rho + d_Btest_rho*betahat_rho;
%%% theta
d_B_theta=d_polybasis(X,poly_d,4);
d_Btest_theta=d_polybasis(Xtest,poly_d,4);
dfB_theta=d_B_theta*betahat_theta;
dYnew_theta=theta - dfB_theta;
theta_SK = w_theta'*dYnew_theta + d_Btest_theta*betahat_theta;




f_enhanced=w_enhanced'*[Ynew;dYnew_delta;dYnew_vega;dYnew_rho;dYnew_theta]+Btest*betahat;
df_delta_enhanced=w_delta_enhanced'*[Ynew;dYnew_delta;dYnew_vega;dYnew_rho;dYnew_theta]+d_Btest_delta*betahat;
df_vega_enhanced=w_vega_enhanced'*[Ynew;dYnew_delta;dYnew_vega;dYnew_rho;dYnew_theta]+d_Btest_vega*betahat;
df_rho_enhanced=w_rho_enhanced'*[Ynew;dYnew_delta;dYnew_vega;dYnew_rho;dYnew_theta]+d_Btest_rho*betahat;
df_theta_enhanced=w_theta_enhanced'*[Ynew;dYnew_delta;dYnew_vega;dYnew_rho;dYnew_theta]+d_Btest_theta*betahat;


RMSE_V.SK = sqrt(mean((f_SK - Y_true).^2));
RMSE_V.ESK = sqrt(mean((f_enhanced - Y_true).^2));
Bias_V.SK = (mean((f_SK - Y_true)));
Bias_V.ESK = (mean((f_enhanced - Y_true)));
std_V.SK = std(f_SK - Y_true);
std_V.ESK = std(f_enhanced - Y_true);



RMSE_delta.SK = sqrt(mean((delta_SK - delta_true).^2));
RMSE_delta.ESK = sqrt(mean((df_delta_enhanced - delta_true).^2));
Bias_delta.SK = (mean((delta_SK - delta_true)));
Bias_delta.ESK = (mean((df_delta_enhanced - delta_true)));
std_delta.SK = std(delta_SK - delta_true);
std_delta.ESK = std(df_delta_enhanced - delta_true);


RMSE_vega.SK = sqrt(mean((vega_SK - vega_true).^2));
RMSE_vega.ESK = sqrt(mean((df_vega_enhanced - vega_true).^2));
Bias_vega.SK = (mean((vega_SK - vega_true)));
Bias_vega.ESK = (mean((df_vega_enhanced - vega_true)));
std_vega.SK = std(vega_SK - vega_true);
std_vega.ESK = std(df_vega_enhanced - vega_true);


RMSE_rho.SK = sqrt(mean((rho_SK - rho_true).^2));
RMSE_rho.ESK = sqrt(mean((df_rho_enhanced - rho_true).^2));
Bias_rho.SK = (mean((rho_SK - rho_true)));
Bias_rho.ESK = (mean((df_rho_enhanced - rho_true)));
std_rho.SK = std(rho_SK - rho_true);
std_rho.ESK = std(df_rho_enhanced - rho_true);


RMSE_theta.SK = sqrt(mean((theta_SK - theta_true).^2));
RMSE_theta.ESK = sqrt(mean((df_theta_enhanced - theta_true).^2));
Bias_theta.SK = (mean((theta_SK - theta_true)));
Bias_theta.ESK = (mean((df_theta_enhanced - theta_true)));
std_theta.SK = std(theta_SK - theta_true);
std_theta.ESK = std(df_theta_enhanced - theta_true);



% subplot(1,2,1)
% plot(Xtest(:,1),Y_true,'r*',Xtest(:,1),f_SK,'ko',Xtest(:,1),f_enhanced,'bd',Xtest(:,1),f_enhanced_const,'g^')
% legend('True','SK','ESK','ESK-const')
% subplot(1,2,2)
% plot(Xtest(:,3),rho_true,'r*',Xtest(:,3),rho_SK,'ko',Xtest(:,3),df_rho_enhanced,'bd',Xtest(:,3),df_rho_enhanced_const,'g^')
% legend('True','SK','ESK','ESK-const')
1
%subplot(1,2,1)

% idx=1;
% dw = SKpredict_M_linear_dw(SKmodel1,Xtest,Btest,idx);
% 
% dw_enhanced = SKpredict_grad1_M_linear_dw(SKmodel_enhanced,Xtest,Btest,idx);
% 
% df2=dw'*Ynew + d_Btest*betahat;
% 
% %df_enhanced = dw_enhanced' * [Ynew;dYnew] + d_Btest*betahat;
% 
% subplot(1,2,1)
% %plot(X,Y,'kx',Xtest,Y_true,'r',Xtest,f1,'b--',Xtest,f2,'c:','LineWidth',2)
% %legend('point','True','Barry','Linearcomb')
% plot(X,Y,'kx',Xtest,Y_true,'r',Xtest,f2,'c--',Xtest,f_enhanced,'b--','LineWidth',2)
% legend('point','True','SK','ESK')
% subplot(1,2,2)
% plot(X,dY,'kx',Xtest,dY_true,'r',Xtest,df_only,'c--',Xtest,df_enhanced,'b--','LineWidth',2)
% legend('point','True','D-SK','ESK-D')


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
    
    

