figure(1)
hold on
boxplot([RMSE_V_SK',RMSE_V_ESK'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],RMSEnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],RMSEaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],RMSEadddc*ones(1,2),'c-')
ylabel('RMSE for Price')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(2)
hold on
boxplot([RMSE_delta_SK',RMSE_delta_ESK'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE for Delta')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(3)
hold on
boxplot([RMSE_vega_SK',RMSE_vega_ESK'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE for Vega')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(4)
hold on
boxplot([RMSE_rho_SK',RMSE_rho_ESK'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE for Rho')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(5)
hold on
boxplot([RMSE_theta_SK',RMSE_theta_ESK'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE for Theta')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(6)
hold on
boxplot([Bias_V_SK',Bias_V_ESK'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],RMSEnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],RMSEaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],RMSEadddc*ones(1,2),'c-')
ylabel('Bias for Price')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(7)
hold on
boxplot([Bias_delta_SK',Bias_delta_ESK'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('Bias for Delta')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(8)
hold on
boxplot([Bias_vega_SK',Bias_vega_ESK'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('Bias for Vega')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(9)
hold on
boxplot([Bias_rho_SK',Bias_rho_ESK'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('Bias for Rho')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(10)
hold on
boxplot([Bias_theta_SK',Bias_theta_ESK'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('Bias for Theta')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(11)
hold on
boxplot([std_V_SK',std_V_ESK'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],RMSEnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],RMSEaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],RMSEadddc*ones(1,2),'c-')
ylabel('Standard Deviation for Price')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(12)
hold on
boxplot([std_delta_SK',std_delta_ESK'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('Standard Deviation for Delta')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(13)
hold on
boxplot([std_vega_SK',std_vega_ESK'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('Standard Deviation for Vega')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(14)
hold on
boxplot([std_rho_SK',std_rho_ESK'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('Standard Deviation for Rho')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(15)
hold on
boxplot([std_theta_SK',std_theta_ESK'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('Standard Deviation for Theta')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off