function [V_c,dV_c]=bls_Euro_call(S,K,r,sigma,t,T)
tau=T-t;
d1=1/(sigma*sqrt(tau))*(log(S./K)+(r+sigma^2/2)*tau);
d2=d1-sigma*sqrt(tau);
V_c = normcdf(d1,0,1).*S-normcdf(d2,0,1).*K*exp(-r*tau);
dV_c= normcdf(d1,0,1);
