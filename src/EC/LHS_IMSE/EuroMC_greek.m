function [V,vV,delta,vdelta,vega,vvega,rho,vrho,theta,vtheta,vmatrix]=EuroMC_greek(S0,r,sigma,T,K,M)
l=length(S0);
V=zeros(l,1);vV=zeros(l,1);delta=zeros(l,1);vdelta=zeros(l,1);
vega=zeros(l,1);vvega=zeros(l,1);rho=zeros(l,1);vrho=zeros(l,1);theta=zeros(l,1);vtheta=zeros(l,1);
vVdelta=zeros(l,1);vVvega=zeros(l,1);vVrho=zeros(l,1);vVT=zeros(l,1);vdeltavega=zeros(l,1);vdeltarho=zeros(l,1);
vdeltatheta=zeros(l,1);vvegarho=zeros(l,1);vvegatheta=zeros(l,1);vrhotheta=zeros(l,1);
for i=1:l
Z=normrnd(0,1,M,1);
ST=S0(i)*exp((r(i)-0.5*sigma(i)^2)*T(i)+sigma(i)*sqrt(T(i))*Z);
Payoff=exp(-r(i)*T(i))*max(ST-K,0);
V(i)=mean(Payoff);
vV(i)=var(Payoff)/M;
%%%delta
IPA_S0=exp(-r(i)*T(i))*ST/S0(i).*(ST>K);
delta(i)=mean(IPA_S0);
vdelta(i)=var(IPA_S0)/M;
%%%vega
IPA_sigma=exp(-r(i)*T(i))*ST.*(-sigma(i)*T(i)+sqrt(T(i))*Z).*(ST>K);
vega(i)=mean(IPA_sigma);
vvega(i)=var(IPA_sigma)/M;
%%%rho
IPA_rho=-T(i)*exp(-r(i)*T(i))*(ST - K).*(ST>K)+ exp(-r(i)*T(i))*T(i)*ST.*(ST>K);
rho(i)=mean(IPA_rho);
vrho(i)=var(IPA_rho)/M;
%%%theta
IPA_T=(-r(i)*exp(-r(i)*T(i))*max(ST-K,0)+exp(-r(i)*T(i))*((ST>K).*(((r(i)-0.5*sigma(i)^2)+0.5*T(i)^(-0.5)*sigma(i)*Z).*ST)));
theta(i)=mean(IPA_T);
vtheta(i)=var(IPA_T)/M;
%%% vmatrix
covVdelta = cov(IPA_S0,Payoff); vVdelta(i) = covVdelta(1,2)/M;
covVvega = cov(IPA_sigma,Payoff);vVvega(i) = covVvega(1,2)/M;
covVrho = cov(IPA_rho,Payoff);vVrho(i) = covVrho(1,2)/M;
covVT = cov(IPA_T,Payoff);vVT(i) = covVT(1,2)/M;
covdeltavega = cov(IPA_S0,IPA_sigma);vdeltavega(i) = covdeltavega(1,2)/M;
covdeltarho = cov(IPA_S0,IPA_rho);vdeltarho(i) = covdeltarho(1,2)/M;
covdeltatheta = cov(IPA_S0, IPA_T);vdeltatheta(i) = covdeltatheta(1,2)/M;
covvegarho = cov(IPA_sigma,IPA_rho);vvegarho(i)=covvegarho(1,2)/M;
covvegatheta = cov(IPA_sigma,IPA_T);vvegatheta(i)=covvegatheta(1,2)/M;
covrhotheta = cov(IPA_rho,IPA_T);vrhotheta(i)=covrhotheta(1,2)/M;
end
vmatrix=[diag(vV),     diag(vVdelta),    diag(vVvega),    diag(vVrho),    diag(vVT);
         diag(vVdelta),diag(vdelta),     diag(vdeltavega),diag(vdeltarho),diag(vdeltatheta);
         diag(vVvega), diag(vdeltavega), diag(vvega),     diag(vvegarho), diag(vvegatheta);
         diag(vVrho),  diag(vdeltarho),  diag(vvegarho),  diag(vrho),     diag(vrhotheta);
         diag(vVT),diag(vdeltatheta),diag(vvegatheta),diag(vrhotheta),diag(vtheta)];


