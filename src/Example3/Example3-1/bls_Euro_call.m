function [V_c,delta_c,sigma_c,theta_c]=bls_Euro_call(S,T,r,sigma,K)
tau=T;
d1=1./(sigma*sqrt(tau)).*(log(S./K)+(r+sigma^2/2)*tau);
d2=d1-sigma*sqrt(tau);
V_c = normcdf(d1,0,1).*S-normcdf(d2,0,1).*K.*exp(-r*tau);
delta_c = normcdf(d1,0,1);
sigma_c = S.*normpdf(d1,0,1).*sqrt(tau);
theta_c = -S.*normpdf(d1,0,1)*sigma./(2*sqrt(tau)) - r*K*exp(-r*tau).*normcdf(d2,0,1);

