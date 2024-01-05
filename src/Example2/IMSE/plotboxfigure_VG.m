figure(1)
hold on
boxplot(sqrt([RMSE_PV_SK',RMSE_PV_ESK']),'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],RMSEnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],RMSEaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],RMSEadddc*ones(1,2),'c-')
ylabel('RMSE for Portfolio price P')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(2)
hold on
boxplot([RMSE_PdS0_SK1',RMSE_PdS0_ESK1'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE for dP/dS_01')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(3)
hold on
boxplot([RMSE_PdS0_SK2',RMSE_PdS0_ESK2'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE for dP/dS_02')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(4)
hold on
boxplot([RMSE_PdS0_SK3',RMSE_PdS0_ESK3'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE for dP/dS_01')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(5)
hold on
boxplot([RMSE_PdS0_SK4',RMSE_PdS0_ESK4'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE for dP/dS_04')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(6)
hold on
boxplot([RMSE_PdS0_SK5',RMSE_PdS0_ESK5'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE for dP/dS_05')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(7)
hold on
boxplot([RMSE_Pdsigma_SK1',RMSE_Pdsigma_ESK1'],'Labels',{'SK','GESK'})
ylabel('RMSE for dP/d\sigma_1')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off


figure(8)
hold on
boxplot([RMSE_Pdsigma_SK2',RMSE_Pdsigma_ESK2'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE for dP/d\sigma_2')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(9)
hold on
boxplot([RMSE_Pdsigma_SK3',RMSE_Pdsigma_ESK3'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE for dP/d\sigma_3')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off


figure(10)
hold on
boxplot([RMSE_Pdsigma_SK4',RMSE_Pdsigma_ESK4'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE for dP/d\sigma_4')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off


figure(11)
hold on
boxplot([RMSE_Pdsigma_SK5',RMSE_Pdsigma_ESK5'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE for dP/d\sigma_5')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(12)
hold on
boxplot([RMSE_Pdtheta_SK1',RMSE_Pdtheta_ESK1'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE for dP/d\theta_1')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(13)
hold on
boxplot([RMSE_Pdtheta_SK2',RMSE_Pdtheta_ESK2'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE for dP/d\theta_2')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(14)
hold on
boxplot([RMSE_Pdtheta_SK3',RMSE_Pdtheta_ESK3'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE for dP/d\theta_3')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(15)
hold on
boxplot([RMSE_Pdtheta_SK4',RMSE_Pdtheta_ESK4'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE for dP/d\theta_4')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(16)
hold on
boxplot([RMSE_Pdtheta_SK5',RMSE_Pdtheta_ESK5'],'Labels',{'SK','GESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE for dP/d\theta_5')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off






