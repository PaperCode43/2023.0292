%% data: 5 stocks APPLE, Facebook, Netflix. Alibaba, Tesla
function [RMSE_PV,RMSE_PdS0,RMSE_Pdsigma,RMSE_Pdtheta]=RMSE_portfolio(Ntrue)
addpath ./datafiles
paraset;
% sigma_v=[0.2636 0.2625 0.4012 0.2842 0.4660];
% nu_v = [0.0387 0.0355 0.0394 0.0017 0.0933];
% theta_v = [-0.5185 -0.8288 -1.2344 -2.6984 -1.1459];
% S0_v=[204.47 144.96 303.47 144.85 350.51];
% r_n=0.0273;
% K_v = [204,140,290,150,360];% one at-the-money, two in the money, two out of money 
% TA=1;TL=1;nA=12;nL=12;
%% Training data
X_AAPL=AAPL_X(S0_v(1),sigma_v(1),theta_v(1));
X_FB=Facebook_X(S0_v(2),sigma_v(2),theta_v(2));
X_NFL=Netflix_X(S0_v(3),sigma_v(3),theta_v(3));
X_Ali=Ali_X(S0_v(4),sigma_v(4),theta_v(4));
X_Tesla=Tesla_X(S0_v(5),sigma_v(5),theta_v(5));

%% Testing data and true value
% AAPL Asian
[Xtest_AAPL_Asian,Ytest_true_AAPL_Asian]=Xtest_Asian(S0_v(1),sigma_v(1),theta_v(1),K_v(1),TA,nA,r_n,nu_v(1),Ntrue);
% AAPL Lookback
[Xtest_AAPL_LB,Ytest_true_AAPL_LB]=Xtest_Lookback(S0_v(1),sigma_v(1),theta_v(1),K_v(1),TL,nL,r_n,nu_v(1),Ntrue);

% Facebook Asian
[Xtest_FB_Asian,Ytest_true_FB_Asian]=Xtest_Asian(S0_v(2),sigma_v(2),theta_v(2),K_v(2),TA,nA,r_n,nu_v(2),Ntrue);
% Facebook Lookback
[Xtest_FB_LB,Ytest_true_FB_LB]=Xtest_Lookback(S0_v(2),sigma_v(2),theta_v(2),K_v(2),TL,nL,r_n,nu_v(2),Ntrue);

% Netflix Asian
[Xtest_NFL_Asian,Ytest_true_NFL_Asian]=Xtest_Asian(S0_v(3),sigma_v(3),theta_v(3),K_v(3),TA,nA,r_n,nu_v(3),Ntrue);
% Netflix Lookback
[Xtest_NFL_LB,Ytest_true_NFL_LB]=Xtest_Lookback(S0_v(3),sigma_v(3),theta_v(3),K_v(3),TL,nL,r_n,nu_v(3),Ntrue);

% Alibaba Asian
[Xtest_Ali_Asian,Ytest_true_Ali_Asian]=Xtest_Asian(S0_v(4),sigma_v(4),theta_v(4),K_v(4),TA,nA,r_n,nu_v(4),Ntrue);
% Alibaba Lookback
[Xtest_Ali_LB,Ytest_true_Ali_LB]=Xtest_Lookback(S0_v(4),sigma_v(4),theta_v(4),K_v(4),TL,nL,r_n,nu_v(4),Ntrue);

% Tesla Asian
[Xtest_Tesla_Asian,Ytest_true_Tesla_Asian]=Xtest_Asian(S0_v(5),sigma_v(5),theta_v(5),K_v(5),TA,nA,r_n,nu_v(5),Ntrue);
% Tesla Lookback
[Xtest_Tesla_LB,Ytest_true_Tesla_LB]=Xtest_Lookback(S0_v(5),sigma_v(5),theta_v(5),K_v(5),TL,nL,r_n,nu_v(5),Ntrue);

save('testdata_P')
%% Learning Procedure
M=500000;poly_d=0;
% AAPL Asian and Lookback
[V_AAPL_Asian,VdS0_AAPL_Asian,Vdsigma_AAPL_Asian,Vdtheta_AAPL_Asian]...
    =KrigEuro_cokriging_Asian_p(X_AAPL,K_v(1),TA,nA,r_n,nu_v(1),M,poly_d,Xtest_AAPL_Asian);
[V_AAPL_LB,VdS0_AAPL_LB,Vdsigma_AAPL_LB,Vdtheta_AAPL_LB]...
    =KrigEuro_cokriging_Lookback_p(X_AAPL,K_v(1),TL,nL,r_n,nu_v(1),M,poly_d,Xtest_AAPL_LB);
% Facebook Asian and Lookback
[V_FB_Asian,VdS0_FB_Asian,Vdsigma_FB_Asian,Vdtheta_FB_Asian]...
    =KrigEuro_cokriging_Asian_p(X_FB,K_v(2),TA,nA,r_n,nu_v(2),M,poly_d,Xtest_FB_Asian);
[V_FB_LB,VdS0_FB_LB,Vdsigma_FB_LB,Vdtheta_FB_LB]...
    =KrigEuro_cokriging_Lookback_p(X_FB,K_v(2),TL,nL,r_n,nu_v(2),M,poly_d,Xtest_FB_LB);
% Netflix Asian and Lookback
[V_NFL_Asian,VdS0_NFL_Asian,Vdsigma_NFL_Asian,Vdtheta_NFL_Asian]...
    =KrigEuro_cokriging_Asian_p(X_NFL,K_v(3),TA,nA,r_n,nu_v(3),M,poly_d,Xtest_NFL_Asian);
[V_NFL_LB,VdS0_NFL_LB,Vdsigma_NFL_LB,Vdtheta_NFL_LB]...
    =KrigEuro_cokriging_Lookback_p(X_NFL,K_v(3),TL,nL,r_n,nu_v(3),M,poly_d,Xtest_NFL_LB);
% Alibaba Asian and Lookback
[V_Ali_Asian,VdS0_Ali_Asian,Vdsigma_Ali_Asian,Vdtheta_Ali_Asian]...
    =KrigEuro_cokriging_Asian_p(X_Ali,K_v(4),TA,nA,r_n,nu_v(4),M,poly_d,Xtest_Ali_Asian);
[V_Ali_LB,VdS0_Ali_LB,Vdsigma_Ali_LB,Vdtheta_Ali_LB]...
    =KrigEuro_cokriging_Lookback_p(X_Ali,K_v(4),TL,nL,r_n,nu_v(4),M,poly_d,Xtest_Ali_LB);
% Tesla Asian and Lookback
[V_Tesla_Asian,VdS0_Tesla_Asian,Vdsigma_Tesla_Asian,Vdtheta_Tesla_Asian]...
    =KrigEuro_cokriging_Asian_p(X_Tesla,K_v(5),TA,nA,r_n,nu_v(5),M,poly_d,Xtest_Tesla_Asian);
[V_Tesla_LB,VdS0_Tesla_LB,Vdsigma_Tesla_LB,Vdtheta_Tesla_LB]...
    =KrigEuro_cokriging_Lookback_p(X_Tesla,K_v(5),TL,nL,r_n,nu_v(5),M,poly_d,Xtest_Tesla_LB);

%% Portfolio: true and estimate
w = ones(10,1); %notice that the returned V is 100x1
% Portfolio value
PV_true = [Ytest_true_AAPL_Asian.V,Ytest_true_AAPL_LB.V,Ytest_true_FB_Asian.V,Ytest_true_FB_LB.V,Ytest_true_NFL_Asian.V,Ytest_true_NFL_LB.V,...
    Ytest_true_Ali_Asian.V,Ytest_true_Ali_LB.V,Ytest_true_Tesla_Asian.V,Ytest_true_Tesla_LB.V] * w;
PV_SK = [V_AAPL_Asian.SK,V_AAPL_LB.SK,V_FB_Asian.SK,V_FB_LB.SK,V_NFL_Asian.SK,V_NFL_LB.SK,...
    V_Ali_Asian.SK,V_Ali_LB.SK,V_Tesla_Asian.SK,V_Tesla_LB.SK] * w;
PV_ESK = [V_AAPL_Asian.ESK,V_AAPL_LB.ESK,V_FB_Asian.ESK,V_FB_LB.ESK,V_NFL_Asian.ESK,V_NFL_LB.ESK,...
    V_Ali_Asian.ESK,V_Ali_LB.ESK,V_Tesla_Asian.ESK,V_Tesla_LB.ESK] * w;
% Portfolio dS0
w2 = ones(2,1);
PdS01_true = [Ytest_true_AAPL_Asian.dS0,Ytest_true_AAPL_LB.dS0] * w2;
PdS02_true = [Ytest_true_FB_Asian.dS0,Ytest_true_FB_LB.dS0] * w2;
PdS03_true = [Ytest_true_NFL_Asian.dS0,Ytest_true_NFL_LB.dS0] * w2;
PdS04_true = [Ytest_true_Ali_Asian.dS0,Ytest_true_Ali_LB.dS0] * w2;
PdS05_true = [Ytest_true_Tesla_Asian.dS0,Ytest_true_Tesla_LB.dS0] * w2;

PdS01_SK = [VdS0_AAPL_Asian.SK,VdS0_AAPL_LB.SK] * w2;
PdS02_SK = [VdS0_FB_Asian.SK,VdS0_FB_LB.SK] * w2;
PdS03_SK = [VdS0_NFL_Asian.SK,VdS0_NFL_LB.SK] * w2;
PdS04_SK = [VdS0_Ali_Asian.SK,VdS0_Ali_LB.SK] * w2;
PdS05_SK = [VdS0_Tesla_Asian.SK,VdS0_Tesla_LB.SK] * w2;

PdS01_ESK = [VdS0_AAPL_Asian.ESK,VdS0_AAPL_LB.ESK] * w2;
PdS02_ESK = [VdS0_FB_Asian.ESK,VdS0_FB_LB.ESK] * w2;
PdS03_ESK = [VdS0_NFL_Asian.ESK,VdS0_NFL_LB.ESK] * w2;
PdS04_ESK = [VdS0_Ali_Asian.ESK,VdS0_Ali_LB.ESK] * w2;
PdS05_ESK = [VdS0_Tesla_Asian.ESK,VdS0_Tesla_LB.ESK] * w2;
% Portfolio dsigma
Pdsigma1_true = [Ytest_true_AAPL_Asian.dsigma,Ytest_true_AAPL_LB.dsigma] * w2;
Pdsigma2_true = [Ytest_true_FB_Asian.dsigma,Ytest_true_FB_LB.dsigma] * w2;
Pdsigma3_true = [Ytest_true_NFL_Asian.dsigma,Ytest_true_NFL_LB.dsigma] * w2;
Pdsigma4_true = [Ytest_true_Ali_Asian.dsigma,Ytest_true_Ali_LB.dsigma] * w2;
Pdsigma5_true = [Ytest_true_Tesla_Asian.dsigma,Ytest_true_Tesla_LB.dsigma] * w2;

Pdsigma1_SK = [Vdsigma_AAPL_Asian.SK,Vdsigma_AAPL_LB.SK] * w2;
Pdsigma2_SK = [Vdsigma_FB_Asian.SK,Vdsigma_FB_LB.SK] * w2;
Pdsigma3_SK = [Vdsigma_NFL_Asian.SK,Vdsigma_NFL_LB.SK] * w2;
Pdsigma4_SK = [Vdsigma_Ali_Asian.SK,Vdsigma_Ali_LB.SK] * w2;
Pdsigma5_SK = [Vdsigma_Tesla_Asian.SK,Vdsigma_Tesla_LB.SK] * w2;

Pdsigma1_ESK = [Vdsigma_AAPL_Asian.ESK,Vdsigma_AAPL_LB.ESK] * w2;
Pdsigma2_ESK = [Vdsigma_FB_Asian.ESK,Vdsigma_FB_LB.ESK] * w2;
Pdsigma3_ESK = [Vdsigma_NFL_Asian.ESK,Vdsigma_NFL_LB.ESK] * w2;
Pdsigma4_ESK = [Vdsigma_Ali_Asian.ESK,Vdsigma_Ali_LB.ESK] * w2;
Pdsigma5_ESK = [Vdsigma_Tesla_Asian.ESK,Vdsigma_Tesla_LB.ESK] * w2;

% Portfolio dtheta
Pdtheta1_true = [Ytest_true_AAPL_Asian.dtheta,Ytest_true_AAPL_LB.dtheta] * w2;
Pdtheta2_true = [Ytest_true_FB_Asian.dtheta,Ytest_true_FB_LB.dtheta] * w2;
Pdtheta3_true = [Ytest_true_NFL_Asian.dtheta,Ytest_true_NFL_LB.dtheta] * w2;
Pdtheta4_true = [Ytest_true_Ali_Asian.dtheta,Ytest_true_Ali_LB.dtheta] * w2;
Pdtheta5_true = [Ytest_true_Tesla_Asian.dtheta,Ytest_true_Tesla_LB.dtheta] * w2;

Pdtheta1_SK = [Vdtheta_AAPL_Asian.SK,Vdtheta_AAPL_LB.SK] * w2;
Pdtheta2_SK = [Vdtheta_FB_Asian.SK,Vdtheta_FB_LB.SK] * w2;
Pdtheta3_SK = [Vdtheta_NFL_Asian.SK,Vdtheta_NFL_LB.SK] * w2;
Pdtheta4_SK = [Vdtheta_Ali_Asian.SK,Vdtheta_Ali_LB.SK] * w2;
Pdtheta5_SK = [Vdtheta_Tesla_Asian.SK,Vdtheta_Tesla_LB.SK] * w2;

Pdtheta1_ESK = [Vdtheta_AAPL_Asian.ESK,Vdtheta_AAPL_LB.ESK] * w2;
Pdtheta2_ESK = [Vdtheta_FB_Asian.ESK,Vdtheta_FB_LB.ESK] * w2;
Pdtheta3_ESK = [Vdtheta_NFL_Asian.ESK,Vdtheta_NFL_LB.ESK] * w2;
Pdtheta4_ESK = [Vdtheta_Ali_Asian.ESK,Vdtheta_Ali_LB.ESK] * w2;
Pdtheta5_ESK = [Vdtheta_Tesla_Asian.ESK,Vdtheta_Tesla_LB.ESK] * w2;


RMSE_PV.SK = sqrt(mean((PV_SK - PV_true).^2));
RMSE_PV.ESK = sqrt(mean((PV_ESK - PV_true).^2));

RMSE_PdS0.SK1 = sqrt(mean((PdS01_SK - PdS01_true).^2));
RMSE_PdS0.ESK1 = sqrt(mean((PdS01_ESK - PdS01_true).^2));
RMSE_PdS0.SK2 = sqrt(mean((PdS02_SK - PdS02_true).^2));
RMSE_PdS0.ESK2 = sqrt(mean((PdS02_ESK - PdS02_true).^2));
RMSE_PdS0.SK3 = sqrt(mean((PdS03_SK - PdS03_true).^2));
RMSE_PdS0.ESK3 = sqrt(mean((PdS03_ESK - PdS03_true).^2));
RMSE_PdS0.SK4 = sqrt(mean((PdS04_SK - PdS04_true).^2));
RMSE_PdS0.ESK4 = sqrt(mean((PdS04_ESK - PdS04_true).^2));
RMSE_PdS0.SK5 = sqrt(mean((PdS05_SK - PdS05_true).^2));
RMSE_PdS0.ESK5 = sqrt(mean((PdS05_ESK - PdS05_true).^2));



RMSE_Pdsigma.SK1 = sqrt(mean((Pdsigma1_SK - Pdsigma1_true).^2));
RMSE_Pdsigma.ESK1 = sqrt(mean((Pdsigma1_ESK - Pdsigma1_true).^2));
RMSE_Pdsigma.SK2 = sqrt(mean((Pdsigma2_SK - Pdsigma2_true).^2));
RMSE_Pdsigma.ESK2 = sqrt(mean((Pdsigma2_ESK - Pdsigma2_true).^2));
RMSE_Pdsigma.SK3 = sqrt(mean((Pdsigma3_SK - Pdsigma3_true).^2));
RMSE_Pdsigma.ESK3 = sqrt(mean((Pdsigma3_ESK - Pdsigma3_true).^2));
RMSE_Pdsigma.SK4 = sqrt(mean((Pdsigma4_SK - Pdsigma4_true).^2));
RMSE_Pdsigma.ESK4 = sqrt(mean((Pdsigma4_ESK - Pdsigma4_true).^2));
RMSE_Pdsigma.SK5 = sqrt(mean((Pdsigma5_SK - Pdsigma5_true).^2));
RMSE_Pdsigma.ESK5 = sqrt(mean((Pdsigma5_ESK - Pdsigma5_true).^2));




RMSE_Pdtheta.SK1 = sqrt(mean((Pdtheta1_SK - Pdtheta1_true).^2));
RMSE_Pdtheta.ESK1 = sqrt(mean((Pdtheta1_ESK - Pdtheta1_true).^2));
RMSE_Pdtheta.SK2 = sqrt(mean((Pdtheta2_SK - Pdtheta2_true).^2));
RMSE_Pdtheta.ESK2 = sqrt(mean((Pdtheta2_ESK - Pdtheta2_true).^2));
RMSE_Pdtheta.SK3 = sqrt(mean((Pdtheta3_SK - Pdtheta3_true).^2));
RMSE_Pdtheta.ESK3 = sqrt(mean((Pdtheta3_ESK - Pdtheta3_true).^2));
RMSE_Pdtheta.SK4 = sqrt(mean((Pdtheta4_SK - Pdtheta4_true).^2));
RMSE_Pdtheta.ESK4 = sqrt(mean((Pdtheta4_ESK - Pdtheta4_true).^2));
RMSE_Pdtheta.SK5 = sqrt(mean((Pdtheta5_SK - Pdtheta5_true).^2));
RMSE_Pdtheta.ESK5 = sqrt(mean((Pdtheta5_ESK - Pdtheta5_true).^2));


