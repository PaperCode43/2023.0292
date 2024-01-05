rep=50;
K=105;
M=10000;
RMSE_V_SK=zeros(1,rep);RMSE_V_ESK=zeros(1,rep);
RMSE_delta_SK=zeros(1,rep);RMSE_delta_ESK=zeros(1,rep);
RMSE_vega_SK=zeros(1,rep);RMSE_vega_ESK=zeros(1,rep);
RMSE_rho_SK=zeros(1,rep);RMSE_rho_ESK=zeros(1,rep);
RMSE_theta_SK=zeros(1,rep);RMSE_theta_ESK=zeros(1,rep);
Bias_V_SK=zeros(1,rep);Bias_V_ESK=zeros(1,rep);
Bias_delta_SK=zeros(1,rep);Bias_delta_ESK=zeros(1,rep);
Bias_vega_SK=zeros(1,rep);Bias_vega_ESK=zeros(1,rep);
Bias_rho_SK=zeros(1,rep);Bias_rho_ESK=zeros(1,rep);
Bias_theta_SK=zeros(1,rep);Bias_theta_ESK=zeros(1,rep);
std_V_SK=zeros(1,rep);std_V_ESK=zeros(1,rep);
std_delta_SK=zeros(1,rep);std_delta_ESK=zeros(1,rep);
std_vega_SK=zeros(1,rep);std_vega_ESK=zeros(1,rep);
std_rho_SK=zeros(1,rep);std_rho_ESK=zeros(1,rep);
std_theta_SK=zeros(1,rep);std_theta_ESK=zeros(1,rep);
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
    [RMSE_V,Bias_V,std_V,RMSE_delta,Bias_delta,std_delta,RMSE_vega,Bias_vega,std_vega,RMSE_rho,Bias_rho,std_rho,RMSE_theta,Bias_theta,std_theta]...
        =KrigEuro_cokriging_all(X,K,M,0,0,0);
    RMSE_V_SK(i) = RMSE_V.SK;RMSE_V_ESK(i) = RMSE_V.ESK;
    RMSE_delta_SK(i) = RMSE_delta.SK;RMSE_delta_ESK(i) = RMSE_delta.ESK;
    RMSE_vega_SK(i) = RMSE_vega.SK;RMSE_vega_ESK(i) = RMSE_vega.ESK;
    RMSE_rho_SK(i) = RMSE_rho.SK;RMSE_rho_ESK(i) = RMSE_rho.ESK;
    RMSE_theta_SK(i) = RMSE_theta.SK;RMSE_theta_ESK(i) = RMSE_theta.ESK;

    Bias_V_SK(i) = Bias_V.SK;Bias_V_ESK(i) = Bias_V.ESK;
    Bias_delta_SK(i) = Bias_delta.SK;Bias_delta_ESK(i) = Bias_delta.ESK;
    Bias_vega_SK(i) = Bias_vega.SK;Bias_vega_ESK(i) = Bias_vega.ESK;
    Bias_rho_SK(i) = Bias_rho.SK;Bias_rho_ESK(i) = Bias_rho.ESK;
    Bias_theta_SK(i) = Bias_theta.SK;Bias_theta_ESK(i) = Bias_theta.ESK;

    std_V_SK(i) = std_V.SK;std_V_ESK(i) = std_V.ESK;
    std_delta_SK(i) = std_delta.SK;std_delta_ESK(i) = std_delta.ESK;
    std_vega_SK(i) = std_vega.SK;std_vega_ESK(i) = std_vega.ESK;
    std_rho_SK(i) = std_rho.SK;std_rho_ESK(i) = std_rho.ESK;
    std_theta_SK(i) = std_theta.SK;std_theta_ESK(i) = std_theta.ESK;
end
save rep_50_all_2_lubound_bias_var_I.mat
plotboxfigure2