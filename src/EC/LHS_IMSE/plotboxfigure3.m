figure(1)
hold on
boxplot(([RMSE_V_SKm',RMSE_V_SKM',RMSE_V_SKml',RMSE_V_SKc',RMSE_V_ESKm',RMSE_V_ESKM',RMSE_V_ESKml',RMSE_V_ESKc']),...
    'Labels',{'SK-mM','SK-Mm','SK-Matlab','SK-Matlab-CR','GESK-mM','GESK-Mm','GESK-Matlab','GESK-Matlab-CR'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],RMSEnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],RMSEaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],RMSEadddc*ones(1,2),'c-')
ylabel('RMSE for Price')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(2)
hold on
boxplot(([RMSE_delta_SKm',RMSE_delta_SKM',RMSE_delta_SKml',RMSE_delta_SKc',RMSE_delta_ESKm',RMSE_delta_ESKM',RMSE_delta_ESKml',RMSE_delta_ESKc']),...
   'Labels',{'SK-mM','SK-Mm','SK-Matlab','SK-Matlab-CR','GESK-mM','GESK-Mm','GESK-Matlab','GESK-Matlab-CR'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE for Delta')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(3)
hold on
boxplot(([RMSE_vega_SKm',RMSE_vega_SKM',RMSE_vega_SKml',RMSE_vega_SKc',RMSE_vega_ESKm',RMSE_vega_ESKM',RMSE_vega_ESKml',RMSE_vega_ESKc']),...
    'Labels',{'SK-mM','SK-Mm','SK-Matlab','SK-Matlab-CR','GESK-mM','GESK-Mm','GESK-Matlab','GESK-Matlab-CR'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE for Vega')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(4)
hold on
boxplot(([RMSE_rho_SKm',RMSE_rho_SKM',RMSE_rho_SKml',RMSE_rho_SKc',RMSE_rho_ESKm',RMSE_rho_ESKM',RMSE_rho_ESKml',RMSE_rho_ESKc']),...
    'Labels',{'SK-mM','SK-Mm','SK-Matlab','SK-Matlab-CR','GESK-mM','GESK-Mm','GESK-Matlab','GESK-Matlab-CR'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE for Rho')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(5)
hold on
boxplot(([RMSE_theta_SKm',RMSE_theta_SKM',RMSE_theta_SKml',RMSE_theta_SKc',RMSE_theta_ESKm',RMSE_theta_ESKM',RMSE_theta_ESKml',RMSE_theta_ESKc']),...
    'Labels',{'SK-mM','SK-Mm','SK-Matlab','SK-Matlab-CR','GESK-mM','GESK-Mm','GESK-Matlab','GESK-Matlab-CR'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE for Theta')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off