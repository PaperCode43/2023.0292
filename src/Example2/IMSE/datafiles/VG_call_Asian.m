%% This file tests VG Asian option
function [VG_c_Asian,c_delta,c_dsigma,c_dtheta,c_dnu] = VG_call_Asian(S_0,K,T,n,N,r_n,sigma_vg,nu_vg,theta_vg)
dt=T/n;
matx=triu(ones(12,12));
matx(1:6,1:6)=0;
%% price
[S_VG,G,Z]=S_VGT(S_0,T,n,N,r_n,sigma_vg,nu_vg,theta_vg);
dis=exp(-r_n*T);
ariave = mean(S_VG(:,7:end),2);
idh = ariave > K;
h = (ariave - K).*idh*dis;
VG_c_Asian = mean(h);
%% dCdS0
IPA_delta = ariave/S_0.*idh*dis;
c_delta=mean(IPA_delta);
%% dCdsigma
dwdsigma = sigma_vg/(1 - theta_vg * nu_vg - sigma_vg^2*nu_vg/2);
dVGdsigma = sqrt(G(:,:)).*Z(:,:);
prefact_sigma = (-dwdsigma*dt+dVGdsigma)*matx;
dariavedsigma = sum(prefact_sigma(:,7:end).*S_VG(:,7:end),2)/6;
IPA_dVdsigma = dis * dariavedsigma.*idh;
c_dsigma = mean(IPA_dVdsigma);
%% dCdtheta
dwdtheta = 1/(1 - theta_vg * nu_vg - sigma_vg^2*nu_vg/2);
dVGdtheta = G(:,:);
prefact_theta = (-dwdtheta*dt+dVGdtheta)*matx;
dariavedtheta = sum(prefact_theta(:,7:end).*S_VG(:,7:end),2)/6;
IPA_dVdtheta = dis * dariavedtheta.*idh;
c_dtheta = mean(IPA_dVdtheta);
% %% dCdnu
% w = -log(1 - theta_vg*nu_vg - sigma_vg^2*nu_vg/2)/nu_vg;
% dwdnu = ((theta_vg+sigma_vg^2/2)/(1 - theta_vg*nu_vg - sigma_vg^2*nu_vg/2) - w)/nu_vg;
% dVGdnu=zeros(N,n);
% for j=1:n
%     j
%     
% Y = G(:,j)/nu_vg;
% %Y = gamrnd(T/nu_vg,1,N,1);
% f1=@(s)s.^(dt/nu_vg - 1).*exp(-s).*log(s);
% I1=integral(f1,0,inf);
% I2 = zeros(N,1);
% %warning('off')
% %%% if dt/nu_vg-1<0, then the integral will be singular at 0, so we need to
% %%% deal with Y
% if dt/nu_vg - 1<0
%     Y=Y+eps;
% end
% for i=1:N
%     
%     I2(i)=integral(f1,eps,Y(i));
% end
% u = gamcdf(Y,dt/nu_vg,1);
% dYdnu = Y.^(1-dt/nu_vg).*exp(Y)*dt/nu_vg^2.*(-u*I1+I2);
% dGTdnu = Y + nu_vg *dYdnu;
% dVGdnu(:,j) = theta_vg * dGTdnu + sigma_vg/2 * (G(:,j)).^(-1/2).*Z(:,j).*dGTdnu;
% 
% end
% prefact_nu = (-dwdnu*dt+dVGdnu)*matx;
% dariavednu = sum(prefact_nu.*S_VG,2)/n;
% IPA_dVdnu = dis * dariavednu.*idh;
% c_dnu = mean(IPA_dVdnu);
c_dnu=1;