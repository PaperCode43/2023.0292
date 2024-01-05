function [C_tr,dVdS0,dVdr,dVdT,dVdsigma,dVdnu,dVdtheta]=VG_call_price_Greeks(S_0,T,r_n,sigma_vg,nu_vg,theta_vg,K)

C_tr = CallPricingFFT_VG(18,S_0,K,T,r_n,0,sigma_vg,nu_vg,theta_vg);
%%% Greeks %%%
epl=0.00001;
%dV/dS0
C_S0_1=CallPricingFFT_VG(18,S_0+epl,K,T,r_n,0,sigma_vg,nu_vg,theta_vg);
C_S0_2=CallPricingFFT_VG(18,S_0-epl,K,T,r_n,0,sigma_vg,nu_vg,theta_vg);
dVdS0=(C_S0_1 - C_S0_2)/(2*epl);
%dV/dr
C_r_1=CallPricingFFT_VG(18,S_0,K,T,r_n+epl,0,sigma_vg,nu_vg,theta_vg);
C_r_2=CallPricingFFT_VG(18,S_0,K,T,r_n-epl,0,sigma_vg,nu_vg,theta_vg);
dVdr=(C_r_1 - C_r_2)/(2*epl);
%dV/dT
C_T_1=CallPricingFFT_VG(18,S_0,K,T+epl,r_n,0,sigma_vg,nu_vg,theta_vg);
C_T_2=CallPricingFFT_VG(18,S_0,K,T-epl,r_n,0,sigma_vg,nu_vg,theta_vg);
dVdT=(C_T_1 - C_T_2)/(2*epl);
%dV/dsigma
C_sigma_1=CallPricingFFT_VG(18,S_0,K,T,r_n,0,sigma_vg+epl,nu_vg,theta_vg);
C_sigma_2=CallPricingFFT_VG(18,S_0,K,T,r_n,0,sigma_vg-epl,nu_vg,theta_vg);
dVdsigma=(C_sigma_1 - C_sigma_2)/(2*epl);
%dV/dnu
C_nu_1=CallPricingFFT_VG(18,S_0,K,T,r_n,0,sigma_vg,nu_vg+epl,theta_vg);
C_nu_2=CallPricingFFT_VG(18,S_0,K,T,r_n,0,sigma_vg,nu_vg-epl,theta_vg);
dVdnu=(C_nu_1 - C_nu_2)/(2*epl);
%dV/dtheta
C_theta_1=CallPricingFFT_VG(18,S_0,K,T,r_n,0,sigma_vg,nu_vg,theta_vg+epl);
C_theta_2=CallPricingFFT_VG(18,S_0,K,T,r_n,0,sigma_vg,nu_vg,theta_vg-epl);
dVdtheta=(C_theta_1 - C_theta_2)/(2*epl);