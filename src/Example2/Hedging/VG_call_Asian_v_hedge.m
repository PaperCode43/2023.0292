%% VG Asian option, part average, the last 3 month 190day to 252day
function [V,vV,dS0,vdS0,vVdS0] = VG_call_Asian_v_hedge(S_0,K,T,iday,n,N,r_n,sigma_vg,nu_vg,theta_vg)
l=length(S_0);
n=12;
T=1;
dt=T/(21*n);
V=zeros(l,1);vV=zeros(l,1);dS0=zeros(l,1);vdS0=zeros(l,1);vVdS0=zeros(l,1);
dis=exp(-r_n*T);
iweek=ceil(iday/21)-1;
for i=1:l
%% price
[S_VG,~,~]=S_VGT(S_0(i),T-dt*iday,(n-iweek),N,r_n,sigma_vg,nu_vg,theta_vg);
ariave = mean(S_VG(:,n-2-iweek:end),2);
idh = ariave > K;
h = (ariave - K).*idh*dis;
V(i) = mean(h);
vV(i) = var(h)/N;
%% dCdS0
IPA_delta = ariave/S_0(i).*idh*dis;
dS0(i)=mean(IPA_delta);
vdS0(i) = var(IPA_delta)/N;
covVdS0 = cov(IPA_delta,h); vVdS0(i) = covVdS0(1,2)/N;
%% Covariance
covVdS0 = cov(IPA_delta,h); vVdS0(i) = covVdS0(1,2)/N;

end


