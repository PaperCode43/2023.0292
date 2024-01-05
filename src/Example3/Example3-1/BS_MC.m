function [Phi,VPhi,delta,Vdelta,Vmatrix]=BS_MC(S0,r,sigma,T,N,K1)
M=size(S0,1);
Phi=zeros(M,1);delta=zeros(M,1);
VPhi=zeros(M,1);Vdelta=zeros(M,1);
vVdS=zeros(M,1);
for i=1:M
    Z0=normrnd(0,1,1,N);
    Z=sigma*Z0;
    ST=S0(i)*exp((r-1/2*sigma^2)*T(i)+sqrt(T(i))*Z);
    p=exp(-r*T(i))*(max(ST-K1,0));
    Phi(i)=mean(p);
    VPhi(i)=var(p)/N;
    
    %%delta
    D=exp(-r*T(i))*(ST>K1).*ST/S0(i);
    delta(i)=mean(D);
    Vdelta(i)=var(D)/N;
%     %%vega
%     v1_MCm=exp(-r*T(i))*ST.*(-sigma*T(i)+sqrt(T(i))*Z0).*(ST>K1);
%     v2_MCm=exp(-r*T(i))*ST.*(-sigma*T(i)+sqrt(T(i))*Z0).*(ST>K2);
%     v3_MCm=exp(-r*T(i))*ST.*(-sigma*T(i)+sqrt(T(i))*Z0).*(ST>K3);
%     v4_MCm=exp(-r*T(i))*ST.*(-sigma*T(i)+sqrt(T(i))*Z0).*(ST>K4);
%     v5_MCm=exp(-r*T(i))*ST.*(-sigma*T(i)+sqrt(T(i))*Z0).*(ST>K5);
%     v6_MCm=exp(-r*T(i))*ST.*(-sigma*T(i)+sqrt(T(i))*Z0).*(ST>K6);
%     v=-2*v2_MCm+v1_MCm+v3_MCm - 2*v5_MCm+v4_MCm+v6_MCm;
%     vega(i)=mean(v);
%     Vvega(i)=var(v)/N;
%%% covariance
    covVdS=cov(p,D);vVdS(i)=covVdS(1,2)/N;

end
Vmatrix=[diag(VPhi),diag(vVdS);diag(vVdS),diag(Vdelta)];