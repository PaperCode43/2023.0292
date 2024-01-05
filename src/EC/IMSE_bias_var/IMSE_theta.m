% S0 test [80,120]
function [IMSE_ESK_Y,IMSE_ESK_dY,IMSE_SK_Y,IMSE_SK_dY]=IMSE_theta(S0,r,sigma,K,M,t,poly_d,rep)
u=lhsdesign(6,rep);
Y_true=zeros(40,rep);%%81 test points 
dY_true=zeros(40,rep);f_enhanced=zeros(40,rep);df_enhanced=zeros(40,rep);
f_SK=zeros(40,rep);df_SK=zeros(40,rep);
IMSE_ESK_Y=zeros(rep,1);IMSE_ESK_dY=zeros(rep,1);
IMSE_SK_Y=zeros(rep,1);IMSE_SK_dY=zeros(rep,1);
for i=1:rep
    X=2*u(:,i);
[Y_true(:,i),dY_true(:,i),f_enhanced(:,i),df_enhanced(:,i),f_SK(:,i),df_SK(:,i)]...
    =KrigEuro_cokriging_theta(S0,r,sigma,X,K,M,t,poly_d,0);
IMSE_ESK_Y(i) = mean((Y_true(:,i) - f_enhanced(:,i)).^2);
IMSE_ESK_dY(i) = mean((dY_true(:,i) - df_enhanced(:,i)).^2);
IMSE_SK_Y(i) = mean((Y_true(:,i) - f_SK(:,i)).^2);
IMSE_SK_dY(i) = mean((dY_true(:,i) - df_SK(:,i)).^2);
end




mIMSE_ESK_Y = mean(IMSE_ESK_Y);
mIMSE_ESK_dY = mean(IMSE_ESK_dY);
mIMSE_SK_Y = mean(IMSE_SK_Y);
mIMSE_SK_dY = mean(IMSE_SK_dY);

figure(1)
hold on
boxplot(sqrt([IMSE_SK_Y,IMSE_ESK_Y]),'Labels',{'SK','ESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],RMSEnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],RMSEaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],RMSEadddc*ones(1,2),'c-')
ylabel('RMSE')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

figure(2)
hold on
boxplot([IMSE_SK_dY,IMSE_ESK_dY],'Labels',{'SK','ESK'})
%plot([0.5+1/4+0.035,1.5-1/4-0.035],PCCnoadd*ones(1,2),'c-',[1.5+0.035+1/4,2.5-1/4-0.035],PCCaddec*ones(1,2),'c-',...
 %   [2.5+0.035+1/4,3.5-1/4-0.035],PCCadddc*ones(1,2),'c-')
ylabel('RMSE')
set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),'FontSize',16)
hold off

%vIMSE_ESK_Y = var(IMSE_ESK_Y);
%vIMSE_ESK_dY = var(IMSE_ESK_dY);
%vIMSE_SK_Y = var(IMSE_SK_Y);
%vIMSE_SK_dY = var(IMSE_SK_dY);