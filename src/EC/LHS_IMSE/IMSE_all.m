rep=50;
K=105;
M=10000;
RMSE_V_SKm=zeros(1,rep);RMSE_V_ESKm=zeros(1,rep);
RMSE_delta_SKm=zeros(1,rep);RMSE_delta_ESKm=zeros(1,rep);
RMSE_vega_SKm=zeros(1,rep);RMSE_vega_ESKm=zeros(1,rep);
RMSE_rho_SKm=zeros(1,rep);RMSE_rho_ESKm=zeros(1,rep);
RMSE_theta_SKm=zeros(1,rep);RMSE_theta_ESKm=zeros(1,rep);
RMSE_V_SKM=zeros(1,rep);RMSE_V_ESKM=zeros(1,rep);
RMSE_delta_SKM=zeros(1,rep);RMSE_delta_ESKM=zeros(1,rep);
RMSE_vega_SKM=zeros(1,rep);RMSE_vega_ESKM=zeros(1,rep);
RMSE_rho_SKM=zeros(1,rep);RMSE_rho_ESKM=zeros(1,rep);
RMSE_theta_SKM=zeros(1,rep);RMSE_theta_ESKM=zeros(1,rep);
RMSE_V_SKc=zeros(1,rep);RMSE_V_ESKc=zeros(1,rep);
RMSE_delta_SKc=zeros(1,rep);RMSE_delta_ESKc=zeros(1,rep);
RMSE_vega_SKc=zeros(1,rep);RMSE_vega_ESKc=zeros(1,rep);
RMSE_rho_SKc=zeros(1,rep);RMSE_rho_ESKc=zeros(1,rep);
RMSE_theta_SKc=zeros(1,rep);RMSE_theta_ESKc=zeros(1,rep);
RMSE_V_SKml=zeros(1,rep);RMSE_V_ESKml=zeros(1,rep);
RMSE_delta_SKml=zeros(1,rep);RMSE_delta_ESKml=zeros(1,rep);
RMSE_vega_SKml=zeros(1,rep);RMSE_vega_ESKml=zeros(1,rep);
RMSE_rho_SKml=zeros(1,rep);RMSE_rho_ESKml=zeros(1,rep);
RMSE_theta_SKml=zeros(1,rep);RMSE_theta_ESKml=zeros(1,rep);
%% generate LHS designs
U_mini = readmatrix('miniLHD.xlsx');
U_Mm = readmatrix('MmLHD.xlsx');


%%
u_mini=zeros(40,4,rep);
u_Mm=zeros(40,4,rep);
u_matlab=zeros(40,4,rep);
for i=1:rep
    u_mini(:,:,i)=U_mini(i:rep:end,:);
    u_Mm(:,:,i)=U_Mm(i:rep:end,:);
end

parfor i=1:rep
    i
    u_m=u_mini(:,:,i);
    u_M=u_Mm(:,:,i);
    %u_c=u_classLHD(:,:,i);
    %u_ml=u_matlab(:,:,i);
    X_m=zeros(40,4);X_M=zeros(40,4);X_c=zeros(40,4);X_ml=zeros(40,4);
    % minimax_LHD
    X_m(:,1) = (1+(2*u_m(:,1)-1)*0.2)*100;%%S0 [80,120]
    X_m(:,2) = u_m(:,2)*0.3;%%sigma [0,0.3]
    X_m(:,3) = u_m(:,3)*0.1;%%r [0,0.1]
    X_m(:,4) = u_m(:,4)*2;%%T[0,2]
    rng(i);
    [RMSE_Vm,RMSE_deltam,RMSE_vegam,RMSE_rhom,RMSE_thetam]=KrigEuro_cokriging_all(X_m,K,M,0,0,0);
    RMSE_V_SKm(i) = RMSE_Vm.SK;RMSE_V_ESKm(i) = RMSE_Vm.ESK;
    RMSE_delta_SKm(i) = RMSE_deltam.SK;RMSE_delta_ESKm(i) = RMSE_deltam.ESK;
    RMSE_vega_SKm(i) = RMSE_vegam.SK;RMSE_vega_ESKm(i) = RMSE_vegam.ESK;
    RMSE_rho_SKm(i) = RMSE_rhom.SK;RMSE_rho_ESKm(i) = RMSE_rhom.ESK;
    RMSE_theta_SKm(i) = RMSE_thetam.SK;RMSE_theta_ESKm(i) = RMSE_thetam.ESK;
    % Maximin_LHD
    X_M(:,1) = (1+(2*u_M(:,1)-1)*0.2)*100;%%S0 [80,120]
    X_M(:,2) = u_M(:,2)*0.3;%%sigma [0,0.3]
    X_M(:,3) = u_M(:,3)*0.1;%%r [0,0.1]
    X_M(:,4) = u_M(:,4)*2;%%T[0,2]
    rng(i);
    [RMSE_VM,RMSE_deltaM,RMSE_vegaM,RMSE_rhoM,RMSE_thetaM]=KrigEuro_cokriging_all(X_M,K,M,0,0,0);
    RMSE_V_SKM(i) = RMSE_VM.SK;RMSE_V_ESKM(i) = RMSE_VM.ESK;
    RMSE_delta_SKM(i) = RMSE_deltaM.SK;RMSE_delta_ESKM(i) = RMSE_deltaM.ESK;
    RMSE_vega_SKM(i) = RMSE_vegaM.SK;RMSE_vega_ESKM(i) = RMSE_vegaM.ESK;
    RMSE_rho_SKM(i) = RMSE_rhoM.SK;RMSE_rho_ESKM(i) = RMSE_rhoM.ESK;
    RMSE_theta_SKM(i) = RMSE_thetaM.SK;RMSE_theta_ESKM(i) = RMSE_thetaM.ESK;
    % Matlab_LHD
    rng(i);
    u_ml=lhsdesign(40,4);
    X_ml(:,1) = (1+(2*u_ml(:,1)-1)*0.2)*100;%%S0 [80,120]
    X_ml(:,2) = u_ml(:,2)*0.3;%%sigma [0,0.3]
    X_ml(:,3) = u_ml(:,3)*0.1;%%r [0,0.1]
    X_ml(:,4) = u_ml(:,4)*2;%%T[0,2]  
    [RMSE_Vml,RMSE_deltaml,RMSE_vegaml,RMSE_rhoml,RMSE_thetaml]=KrigEuro_cokriging_all(X_ml,K,M,0,0,0);
    RMSE_V_SKml(i) = RMSE_Vml.SK;RMSE_V_ESKml(i) = RMSE_Vml.ESK;
    RMSE_delta_SKml(i) = RMSE_deltaml.SK;RMSE_delta_ESKml(i) = RMSE_deltaml.ESK;
    RMSE_vega_SKml(i) = RMSE_vegaml.SK;RMSE_vega_ESKml(i) = RMSE_vegaml.ESK;
    RMSE_rho_SKml(i) = RMSE_rhoml.SK;RMSE_rho_ESKml(i) = RMSE_rhoml.ESK;
    RMSE_theta_SKml(i) = RMSE_thetaml.SK;RMSE_theta_ESKml(i) = RMSE_thetaml.ESK;
    % Matlab_LHD_correlation
    rng(i);
    u_c=lhsdesign(40,4,'Criterion','correlation');
    X_c(:,1) = (1+(2*u_c(:,1)-1)*0.2)*100;%%S0 [80,120]
    X_c(:,2) = u_c(:,2)*0.3;%%sigma [0,0.3]
    X_c(:,3) = u_c(:,3)*0.1;%%r [0,0.1]
    X_c(:,4) = u_c(:,4)*2;%%T[0,2] 
    [RMSE_Vc,RMSE_deltac,RMSE_vegac,RMSE_rhoc,RMSE_thetac]=KrigEuro_cokriging_all(X_c,K,M,0,0,0);
    RMSE_V_SKc(i) = RMSE_Vc.SK;RMSE_V_ESKc(i) = RMSE_Vc.ESK;
    RMSE_delta_SKc(i) = RMSE_deltac.SK;RMSE_delta_ESKc(i) = RMSE_deltac.ESK;
    RMSE_vega_SKc(i) = RMSE_vegac.SK;RMSE_vega_ESKc(i) = RMSE_vegac.ESK;
    RMSE_rho_SKc(i) = RMSE_rhoc.SK;RMSE_rho_ESKc(i) = RMSE_rhoc.ESK;
    RMSE_theta_SKc(i) = RMSE_thetac.SK;RMSE_theta_ESKc(i) = RMSE_thetac.ESK;
    
end
%save rep_50_LHDcompare5mat.mat
plotboxfigure3