Rep=20;
vRMSE_V_ESK=zeros(Rep,1);vRMSE_V_ESKpde=zeros(Rep,1);
vRMSE_theta_ESK=zeros(Rep,1);vRMSE_theta_ESKpde=zeros(Rep,1);
vRMSE_delta_ESK=zeros(Rep,1);vRMSE_delta_ESKpde=zeros(Rep,1);
vRMSE_gamma_ESK=zeros(Rep,1);vRMSE_gamma_ESKpde=zeros(Rep,1);
for k=1:Rep
    k
   [RMSE_V,RMSE_theta,RMSE_delta,RMSE_gamma]=portfolio_BS_kriging_PDE(1,0,k,32,8);%%
   vRMSE_V_ESK(k)=RMSE_V.ESK;vRMSE_V_ESKpde(k)=RMSE_V.ESKpde;
   vRMSE_theta_ESK(k)=RMSE_theta.ESK;vRMSE_theta_ESKpde(k)=RMSE_theta.ESKpde;
   vRMSE_delta_ESK(k)=RMSE_delta.ESK;vRMSE_delta_ESKpde(k)=RMSE_delta.ESKpde;
   vRMSE_gamma_ESK(k)=RMSE_gamma.ESK;vRMSE_gamma_ESKpde(k)=RMSE_gamma.ESKpde;
end
save('Rep20_BS')

figure(1)
ylabel('RMSE of Value')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
boxplot([vRMSE_V_ESK,vRMSE_V_ESKpde],'Labels',{'GESK','PDE-GESK'})

figure(2)
set(gca,'FontSize',16);
ylabel('RMSE of Delta')
set(findobj(gca,'Type','text'),'FontSize',16)
boxplot([vRMSE_delta_ESK,vRMSE_delta_ESKpde],'Labels',{'GESK','PDE-GESK'})

figure(3)
set(gca,'FontSize',16);
ylabel('RMSE of Gamma')
set(findobj(gca,'Type','text'),'FontSize',16)
boxplot([vRMSE_gamma_ESK,vRMSE_gamma_ESKpde],'Labels',{'GESK','PDE-GESK'})

figure(4)
set(gca,'FontSize',16);
ylabel('RMSE of Theta')
set(findobj(gca,'Type','text'),'FontSize',16)
boxplot([vRMSE_theta_ESK,vRMSE_theta_ESKpde],'Labels',{'GESK','PDE-GESK'})






