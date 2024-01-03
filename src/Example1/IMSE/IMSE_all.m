rep=50;
K=105;
M=10000;
RMSE_V_SK=zeros(1,rep);RMSE_V_ESK=zeros(1,rep);
RMSE_delta_SK=zeros(1,rep);RMSE_delta_ESK=zeros(1,rep);
RMSE_vega_SK=zeros(1,rep);RMSE_vega_ESK=zeros(1,rep);
RMSE_rho_SK=zeros(1,rep);RMSE_rho_ESK=zeros(1,rep);
RMSE_theta_SK=zeros(1,rep);RMSE_theta_ESK=zeros(1,rep);
%U=lhsdesign(20,4*rep);
U=sobolset(4*rep);
%u=p(2:21,:);
parfor i=1:rep
    i
    rng(i);
    %u=U(2:41,(i-1)*4+1:i*4);
    u=lhsdesign(40,4);
    X=zeros(40,4);
    X(:,1) = (1+(2*u(:,1)-1)*0.2)*100;%%S0 [80,120]
    X(:,2) = u(:,2)*0.3;%%sigma [0,0.3]
    X(:,3) = u(:,3)*0.1;%%r [0,0.1]
    X(:,4) = u(:,4)*2;%%T[0,2]
    [RMSE_V,RMSE_delta,RMSE_vega,RMSE_rho,RMSE_theta]=KrigEuro_cokriging_all(X,K,M,0,0,0);
    RMSE_V_SK(i) = RMSE_V.SK;RMSE_V_ESK(i) = RMSE_V.ESK;
    RMSE_delta_SK(i) = RMSE_delta.SK;RMSE_delta_ESK(i) = RMSE_delta.ESK;
    RMSE_vega_SK(i) = RMSE_vega.SK;RMSE_vega_ESK(i) = RMSE_vega.ESK;
    RMSE_rho_SK(i) = RMSE_rho.SK;RMSE_rho_ESK(i) = RMSE_rho.ESK;
    RMSE_theta_SK(i) = RMSE_theta.SK;RMSE_theta_ESK(i) = RMSE_theta.ESK;
end
save rep_50.mat
plotboxfigure