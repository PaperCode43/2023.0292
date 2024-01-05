Rep=20;
vRMSE_V_ESK=zeros(Rep,1);vRMSE_V_ESKpde=zeros(Rep,1);
vRMSE_theta_ESK=zeros(Rep,1);vRMSE_theta_ESKpde=zeros(Rep,1);
vARMSE_delta1_ESK=zeros(Rep,1);vARMSE_delta1_ESKpde=zeros(Rep,1);
vARMSE_delta2_ESK=zeros(Rep,1);vARMSE_delta2_ESKpde=zeros(Rep,1);
vARMSE_gamma11_ESK=zeros(Rep,1);vARMSE_gamma11_ESKpde=zeros(Rep,1);
vARMSE_gamma12_ESK=zeros(Rep,1);vARMSE_gamma12_ESKpde=zeros(Rep,1);
vARMSE_gamma22_ESK=zeros(Rep,1);vARMSE_gamma22_ESKpde=zeros(Rep,1);

vRMSE_delta1_ESK=zeros(40,Rep);vRMSE_delta1_ESKpde=zeros(40,Rep);
vRMSE_delta2_ESK=zeros(40,Rep);vRMSE_delta2_ESKpde=zeros(40,Rep);
vRMSE_gamma11_ESK=zeros(40,Rep);vRMSE_gamma11_ESKpde=zeros(40,Rep);
vRMSE_gamma12_ESK=zeros(40,Rep);vRMSE_gamma12_ESKpde=zeros(40,Rep);
vRMSE_gamma22_ESK=zeros(40,Rep);vRMSE_gamma22_ESKpde=zeros(40,Rep);

for k=1:Rep
    k
   [RMSE_V,RMSE_theta,ARMSE_delta1,ARMSE_delta2,ARMSE_gamma11,ARMSE_gamma12,ARMSE_gamma21,ARMSE_gamma22, RMSE_delta1, RMSE_delta2,RMSE_gamma11,RMSE_gamma12,RMSE_gamma22]...
       =portfolio_basket_kriging_PDE_large40(1,0,k,20,10);%%
   %% 50 designs 
   %[RMSE_V,RMSE_theta,ARMSE_delta1,ARMSE_delta2,ARMSE_gamma11,ARMSE_gamma12,ARMSE_gamma21,ARMSE_gamma22, RMSE_delta1, RMSE_delta2,RMSE_gamma11,RMSE_gamma12,RMSE_gamma22]...
       %=portfolio_basket_kriging_PDE_large40(1,0,k,30,20);%%
   vRMSE_V_ESK(k)=RMSE_V.ESK;vRMSE_V_ESKpde(k)=RMSE_V.ESKpde;
   vRMSE_theta_ESK(k)=RMSE_theta.ESK;vRMSE_theta_ESKpde(k)=RMSE_theta.ESKpde;
   vARMSE_delta1_ESK(k)=ARMSE_delta1.ESK;vARMSE_delta1_ESKpde(k)=ARMSE_delta1.ESKpde;
   vARMSE_delta2_ESK(k)=ARMSE_delta2.ESK;vARMSE_delta2_ESKpde(k)=ARMSE_delta2.ESKpde;
   vARMSE_gamma11_ESK(k)=ARMSE_gamma11.ESK;vARMSE_gamma11_ESKpde(k)=ARMSE_gamma11.ESKpde;
   vARMSE_gamma12_ESK(k)=ARMSE_gamma12.ESK;vARMSE_gamma12_ESKpde(k)=ARMSE_gamma12.ESKpde;
   vARMSE_gamma22_ESK(k)=ARMSE_gamma22.ESK;vARMSE_gamma22_ESKpde(k)=ARMSE_gamma22.ESKpde;
   
   vRMSE_delta1_ESK(:,k)=RMSE_delta1.ESK;vRMSE_delta1_ESKpde(:,k)=RMSE_delta1.ESKpde;
   vRMSE_delta2_ESK(:,k)=RMSE_delta2.ESK;vRMSE_delta2_ESKpde(:,k)=RMSE_delta2.ESKpde;
   vRMSE_gamma11_ESK(:,k)=RMSE_gamma11.ESK;vRMSE_gamma11_ESKpde(:,k)=RMSE_gamma11.ESKpde;
   vRMSE_gamma12_ESK(:,k)=RMSE_gamma12.ESK;vRMSE_gamma12_ESKpde(:,k)=RMSE_gamma12.ESKpde;
   vRMSE_gamma22_ESK(:,k)=RMSE_gamma22.ESK;vRMSE_gamma22_ESKpde(:,k)=RMSE_gamma22.ESKpde;
end
%save('Rep20_Port40')

figure(1)
hold on
ylabel('RMSE of Value')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
boxplot([vRMSE_V_ESK,vRMSE_V_ESKpde],'Labels',{'GESK','PDE-GESK'})
hold off

figure(2)
hold on
set(gca,'FontSize',16);
ylabel('Averaged RMSE of Deltas')
set(findobj(gca,'Type','text'),'FontSize',16)
boxplot([vARMSE_delta1_ESK+vARMSE_delta2_ESK,vARMSE_delta1_ESKpde+vARMSE_delta2_ESKpde],'Labels',{'GESK','PDE-GESK'})
hold off

figure(3)
hold on
set(gca,'FontSize',16);
ylabel('Averaged RMSE of Gammas')
set(findobj(gca,'Type','text'),'FontSize',16)
boxplot([vARMSE_gamma11_ESK+2*vARMSE_gamma12_ESK+vARMSE_gamma22_ESK,vARMSE_gamma11_ESKpde+2*vARMSE_gamma12_ESKpde+vARMSE_gamma22_ESKpde],'Labels',{'GESK','PDE-GESK'})
hold off

figure(4)
hold on
set(gca,'FontSize',16);
ylabel('RMSE of Theta')
set(findobj(gca,'Type','text'),'FontSize',16)
boxplot([vRMSE_theta_ESK,vRMSE_theta_ESKpde],'Labels',{'GESK','PDE-GESK'})
hold off





