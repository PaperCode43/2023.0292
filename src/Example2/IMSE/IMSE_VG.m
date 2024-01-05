
Ntrue=1000000;
rep=50;
RMSE_PV_SK=zeros(1,rep);RMSE_PV_ESK=zeros(1,rep);
RMSE_PdS0_SK1=zeros(1,rep);RMSE_PdS0_ESK1=zeros(1,rep);
RMSE_PdS0_SK2=zeros(1,rep);RMSE_PdS0_ESK2=zeros(1,rep);
RMSE_PdS0_SK3=zeros(1,rep);RMSE_PdS0_ESK3=zeros(1,rep);
RMSE_PdS0_SK4=zeros(1,rep);RMSE_PdS0_ESK4=zeros(1,rep);
RMSE_PdS0_SK5=zeros(1,rep);RMSE_PdS0_ESK5=zeros(1,rep);
RMSE_Pdsigma_SK1=zeros(1,rep);RMSE_Pdsigma_ESK1=zeros(1,rep);
RMSE_Pdsigma_SK2=zeros(1,rep);RMSE_Pdsigma_ESK2=zeros(1,rep);
RMSE_Pdsigma_SK3=zeros(1,rep);RMSE_Pdsigma_ESK3=zeros(1,rep);
RMSE_Pdsigma_SK4=zeros(1,rep);RMSE_Pdsigma_ESK4=zeros(1,rep);
RMSE_Pdsigma_SK5=zeros(1,rep);RMSE_Pdsigma_ESK5=zeros(1,rep);
RMSE_Pdtheta_SK1=zeros(1,rep);RMSE_Pdtheta_ESK1=zeros(1,rep);
RMSE_Pdtheta_SK2=zeros(1,rep);RMSE_Pdtheta_ESK2=zeros(1,rep);
RMSE_Pdtheta_SK3=zeros(1,rep);RMSE_Pdtheta_ESK3=zeros(1,rep);
RMSE_Pdtheta_SK4=zeros(1,rep);RMSE_Pdtheta_ESK4=zeros(1,rep);
RMSE_Pdtheta_SK5=zeros(1,rep);RMSE_Pdtheta_ESK5=zeros(1,rep);
for i=1:rep
    i
    rng(i)
    [RMSE_PV,RMSE_PdS0,RMSE_Pdsigma,RMSE_Pdtheta]=RMSE_portfolio(Ntrue);
    RMSE_PV_SK(i) = RMSE_PV.SK; RMSE_PV_ESK(i)=RMSE_PV.ESK;
    
    RMSE_PdS0_SK1(i) = RMSE_PdS0.SK1; RMSE_PdS0_ESK1(i)=RMSE_PdS0.ESK1;
    RMSE_PdS0_SK2(i) = RMSE_PdS0.SK2; RMSE_PdS0_ESK2(i)=RMSE_PdS0.ESK2;
    RMSE_PdS0_SK3(i) = RMSE_PdS0.SK3; RMSE_PdS0_ESK3(i)=RMSE_PdS0.ESK3;
    RMSE_PdS0_SK4(i) = RMSE_PdS0.SK4; RMSE_PdS0_ESK4(i)=RMSE_PdS0.ESK4;
    RMSE_PdS0_SK5(i) = RMSE_PdS0.SK5; RMSE_PdS0_ESK5(i)=RMSE_PdS0.ESK5;
    
    RMSE_Pdsigma_SK1(i) = RMSE_Pdsigma.SK1; RMSE_Pdsigma_ESK1(i)=RMSE_Pdsigma.ESK1;
    RMSE_Pdsigma_SK2(i) = RMSE_Pdsigma.SK2; RMSE_Pdsigma_ESK2(i)=RMSE_Pdsigma.ESK2;
    RMSE_Pdsigma_SK3(i) = RMSE_Pdsigma.SK3; RMSE_Pdsigma_ESK3(i)=RMSE_Pdsigma.ESK3;
    RMSE_Pdsigma_SK4(i) = RMSE_Pdsigma.SK4; RMSE_Pdsigma_ESK4(i)=RMSE_Pdsigma.ESK4;
    RMSE_Pdsigma_SK5(i) = RMSE_Pdsigma.SK5; RMSE_Pdsigma_ESK5(i)=RMSE_Pdsigma.ESK5;
    
    RMSE_Pdtheta_SK1(i) = RMSE_Pdtheta.SK1; RMSE_Pdtheta_ESK1(i)=RMSE_Pdtheta.ESK1;
    RMSE_Pdtheta_SK2(i) = RMSE_Pdtheta.SK2; RMSE_Pdtheta_ESK2(i)=RMSE_Pdtheta.ESK2;
    RMSE_Pdtheta_SK3(i) = RMSE_Pdtheta.SK3; RMSE_Pdtheta_ESK3(i)=RMSE_Pdtheta.ESK3;
    RMSE_Pdtheta_SK4(i) = RMSE_Pdtheta.SK4; RMSE_Pdtheta_ESK4(i)=RMSE_Pdtheta.ESK4;
    RMSE_Pdtheta_SK5(i) = RMSE_Pdtheta.SK5; RMSE_Pdtheta_ESK5(i)=RMSE_Pdtheta.ESK5;
end
paraset;
%save rep_50_VG_portfolio_30design.mat
plotboxfigure_VG