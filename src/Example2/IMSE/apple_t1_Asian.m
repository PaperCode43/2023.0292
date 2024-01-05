S0=139.14;
sigma_vg=0.22;
theta_vg=-0.16;

T=1;nu_vg=0.4339;M_true=1000000;K=140;r_n=0.0103;n=12;

M=10000;poly_d=0;

X=AAPL_X(S0,sigma_vg,theta_vg);
%[Xtest,Ytest_true]=Xtest_Asian(S0,sigma_vg,theta_vg,K,T,n,r_n,nu_vg,M_true);
%save('tempdata1_aapl_Asian')
[Xtest,Ytest_true]=Xtest_Lookback(S0,sigma_vg,theta_vg,K,T,n,r_n,nu_vg,M_true);
save('tempdata1_aapl_Lookback')
[RMSE_V,RMSE_dS0,RMSE_dsigma,RMSE_dtheta]=KrigEuro_cokriging_Asian(X,K,T,r_n,nu_vg,M,poly_d,Xtest,Ytest_true)