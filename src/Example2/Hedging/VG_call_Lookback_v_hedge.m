%% This file tests VG Asian option
function [V,vV,dS0,vdS0,vVdS0] = VG_call_Lookback_v_hedge(S_0,K,T,iday,n,N,r_n,sigma_vg,nu_vg,theta_vg)
l=length(S_0);
V=zeros(l,1);vV=zeros(l,1);dS0=zeros(l,1);vdS0=zeros(l,1);vVdS0=zeros(l,1);
T=1;
n=12;
dt=T/(21*n);
dis=exp(-r_n*T);
iweek=ceil(iday/21)-1;
for i=1:l
%% price
[S_VG,~,~]=S_VGT(S_0(i),T-dt*iday,(n-iweek),N,r_n,sigma_vg,nu_vg,theta_vg);
Smax=max(S_VG(:,n-2-iweek:end),[],2);
idh = Smax> K;
h = dis * (Smax - K).*idh;
V(i) = mean(h);
vV(i) = var(h)/N;
%% dCdS0
dJdS0=Smax/S_0(i);
IPA_delta = dis*dJdS0.*idh;
dS0(i)=mean(IPA_delta);
vdS0(i) = var(IPA_delta)/N;
covVdS0 = cov(IPA_delta,h); vVdS0(i) = covVdS0(1,2)/N;
end
