%% Basket option example test
function [c,c_d1,c_d2,c_d1d1,c_d2d2,c_d1d2,c_d2d1,c_dT]=basket4(S0,T,r,sigma,K,Rho)
%r=0.02;
%S0=[50,70];
%sigma=[0.15,0.3];
%K=60;
%rho=[0.1,0.3];
%% geometric call formula
sigma_new=1/2*sqrt(sigma*Rho*sigma');
F=prod(S0)^(1/2)*exp((r-1/(2*2)*sigma*sigma'+1/2*sigma_new^2)*T);
%F=40:1:100;
d1=(log(F/K)+1/2*sigma_new^2*T)/(sigma_new*sqrt(T));
d2=d1-sigma_new*sqrt(T);
c=exp(-r*T)*(F*normcdf(d1) - K*normcdf(d2));

%%delta
Fd1=1/2*F/S0(1);
c_d1=exp(-r*T)*(Fd1.*normcdf(d1) + F.*normpdf(d1)*(1/(sigma_new*sqrt(T))*1./F.*Fd1) - K*normpdf(d2)*(1/(sigma_new*sqrt(T))*1./F.*Fd1));
Fd2=1/2*F/S0(2);
c_d2=exp(-r*T)*(Fd2.*normcdf(d1) + F.*normpdf(d1)*(1/(sigma_new*sqrt(T))*1./F.*Fd2) - K*normpdf(d2)*(1/(sigma_new*sqrt(T))*1./F.*Fd2));
%%gamma
Fd1d1=1/2*Fd1*S0(1)^(-1)-1/2*F*S0(1)^(-2);
d1dS1=1/(sigma_new*sqrt(T))*1./F.*Fd1;
d1dS1dS1=-1/2*1/(sigma_new*sqrt(T))*S0(1)^(-2);
c_d1d1=exp(-r*T)*( Fd1d1*normcdf(d1) + Fd1*normpdf(d1)*d1dS1 - K*(-d2)*normpdf(d2)*d1dS1^2 - K*normpdf(d2)*d1dS1dS1 + Fd1*normpdf(d1)*d1dS1 + F*(-d1)*normpdf(d1)*d1dS1^2 + F*normpdf(d1)* d1dS1dS1);

Fd2d2=1/2*Fd2*S0(2)^(-1)-1/2*F*S0(2)^(-2);
d1dS2=1/(sigma_new*sqrt(T))*1./F.*Fd2;
d1dS2dS2=-1/2*1/(sigma_new*sqrt(T))*S0(2)^(-2);
c_d2d2=exp(-r*T)*( Fd2d2*normcdf(d1) + Fd2*normpdf(d1)*d1dS2 - K*(-d2)*normpdf(d2)*d1dS2^2 - K*normpdf(d2)*d1dS2dS2 + Fd2*normpdf(d1)*d1dS2 + F*(-d1)*normpdf(d1)*d1dS2^2 + F*normpdf(d1)* d1dS2dS2);

Fd1d2=1/2*Fd2*S0(1)^(-1);
c_d1d2=exp(-r*T)*( Fd1d2*normcdf(d1) + Fd1*normpdf(d1)*d1dS2 - K*(-d2)*normpdf(d2)*d1dS1*d1dS2 + Fd2*normpdf(d1)*d1dS1 + F*(-d1)*normpdf(d1)*d1dS1*d1dS2);

Fd2d1=1/2*Fd1*S0(2)^(-1);
c_d2d1=exp(-r*T)*( Fd2d1*normcdf(d1) + Fd2*normpdf(d1)*d1dS1 - K*(-d2)*normpdf(d2)*d1dS2*d1dS1 + Fd1*normpdf(d1)*d1dS2 + F*(-d1)*normpdf(d1)*d1dS2*d1dS1);

%%Theta
d1dT=(1/2*sigma_new^3*sqrt(T) - (log(F/K)+1/2*sigma_new^2*T)*(sigma_new*1/2*T^(-1/2)))/(sigma_new^2*T);
d2dT=d1dT - 1/2*sigma_new*T^(-1/2);
dFdT=(r-1/(2*2)*sigma*sigma'+1/2*sigma_new^2)*F;
c_dT = -r*c+exp(-r*T)*(F*normpdf(d1)*d1dT + dFdT * normcdf(d1)  - K*normpdf(d2)*d2dT);






