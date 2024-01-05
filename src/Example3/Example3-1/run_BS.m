Rep=20;
vRMSE_V_ESK=zeros(Rep,1);vRMSE_V_ESKpde=zeros(Rep,1);
vRMSE_delta_ESK=zeros(Rep,1);vRMSE_delta_ESKpde=zeros(Rep,1);
for k=1:Rep
    k
   [RMSE_V,RMSE_delta]=portfolio_BS_kriging_PDE(1,0,k,15,5);%%
   vRMSE_V_ESK(k)=RMSE_V.ESK;vRMSE_V_ESKpde(k)=RMSE_V.ESKpde;  
   vRMSE_delta_ESK(k)=RMSE_delta.ESK;vRMSE_delta_ESKpde(k)=RMSE_delta.ESKpde;
end
%save('Rep20_BS')

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







