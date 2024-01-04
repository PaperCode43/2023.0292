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
boxplot([vRMSE_delta_ESK,vRMSE_delta_ESKpde],'Labels',{'GESK','PDE-GESK'})
hold off

figure(3)
hold on
set(gca,'FontSize',16);
ylabel('Averaged RMSE of Gammas')
set(findobj(gca,'Type','text'),'FontSize',16)
boxplot([vRMSE_gamma_ESK,vRMSE_gamma_ESKpde],'Labels',{'GESK','PDE-GESK'})
hold off

figure(4)
hold on
set(gca,'FontSize',16);
ylabel('RMSE of Theta')
set(findobj(gca,'Type','text'),'FontSize',16)
boxplot([vRMSE_theta_ESK,vRMSE_theta_ESKpde],'Labels',{'GESK','PDE-GESK'})
hold off