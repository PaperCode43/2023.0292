%% This file tests VG Asian option Part average 7:12
function [V,vV,dS0,vdS0,dsigma,vdsigma,dtheta,vdtheta,vmatrix] = VG_call_Asian_v(S_0,K,T,n,N,r_n,sigma_vg,nu_vg,theta_vg)
l=length(S_0);
V=zeros(l,1);vV=zeros(l,1);dS0=zeros(l,1);vdS0=zeros(l,1);
dsigma=zeros(l,1);vdsigma=zeros(l,1);dtheta=zeros(l,1);vdtheta=zeros(l,1);
vVdS0=zeros(l,1);vVdsigma=zeros(l,1);vVdtheta=zeros(l,1);vdS0dsigma=zeros(l,1);vdS0dtheta=zeros(l,1);
vdsigmadtheta=zeros(l,1);
dt=T/n;
matx=triu(ones(12,12));
matx(1:6,1:6)=0;
dis=exp(-r_n*T);
for i=1:l
%% price
[S_VG,G,Z]=S_VGT(S_0(i),T,n,N,r_n,sigma_vg(i),nu_vg,theta_vg(i));
ariave = mean(S_VG(:,7:end),2);
idh = ariave > K;
h = (ariave - K).*idh*dis;
V(i) = mean(h);
vV(i) = var(h)/N;
%% dCdS0
IPA_delta = ariave/S_0(i).*idh*dis;
dS0(i)=mean(IPA_delta);
vdS0(i) = var(IPA_delta)/N;
%% dCdsigma
dwdsigma = sigma_vg(i)/(1 - theta_vg(i) * nu_vg - sigma_vg(i)^2*nu_vg/2);
dVGdsigma = sqrt(G(:,:)).*Z(:,:);
prefact_sigma = (-dwdsigma*dt+dVGdsigma)*matx;
dariavedsigma = sum(prefact_sigma(:,7:end).*S_VG(:,7:end),2)/6;
IPA_dVdsigma = dis * dariavedsigma.*idh;
dsigma(i) = mean(IPA_dVdsigma);
vdsigma(i) = var(IPA_dVdsigma)/N;
%% dCdtheta
dwdtheta = 1/(1 - theta_vg(i) * nu_vg - sigma_vg(i)^2*nu_vg/2);
dVGdtheta = G(:,:);
prefact_theta = (-dwdtheta*dt+dVGdtheta)*matx;
dariavedtheta = sum(prefact_theta(:,7:end).*S_VG(:,7:end),2)/6;
IPA_dVdtheta = dis * dariavedtheta.*idh;
dtheta(i) = mean(IPA_dVdtheta);
vdtheta(i) = var(IPA_dVdtheta)/N;
%% Covariance
covVdS0 = cov(IPA_delta,h); vVdS0(i) = covVdS0(1,2)/N;
covVdsigma = cov(IPA_dVdsigma,h); vVdsigma(i) = covVdsigma(1,2)/N;
covVdtheta = cov(IPA_dVdtheta,h); vVdtheta(i) = covVdtheta(1,2)/N;
covdS0dsigma = cov(IPA_delta,IPA_dVdsigma); vdS0dsigma(i) = covdS0dsigma(1,2)/N;
covdS0dtheta = cov(IPA_delta,IPA_dVdtheta); vdS0dtheta(i) = covdS0dtheta(1,2)/N;
covdsigmadtheta = cov(IPA_dVdsigma,IPA_dVdtheta); vdsigmadtheta(i) = covdsigmadtheta(1,2)/N;
end
vmatrix=[diag(vV),       diag(vVdS0),      diag(vVdsigma),       diag(vVdtheta)
         diag(vVdS0),    diag(vdS0),       diag(vdS0dsigma),     diag(vdS0dtheta);
         diag(vVdsigma), diag(vdS0dsigma), diag(vdsigma),        diag(vdsigmadtheta);
         diag(vVdtheta),  diag(vdS0dtheta),diag(vdsigmadtheta),  diag(vdtheta)       ];
