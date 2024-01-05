%% this file to calculate the analytical formula of portfolio value
function [Phi_tau,Phi_tau_delta,Phi_tau_gamma,Phi_tau_theta]=portfolio_basket_price2(S_tau,r,sigma,T,Rho,K1,K2,K3)
%Rho=rho'*rho;
%Rho(logical(eye(size(Rho))))=1;
%sigma_new=1/2*sqrt(sigma*Rho*sigma');
Nt=size(S_tau,1);
Phi_tau=zeros(Nt,1);
Phi_tau_delta=zeros(Nt,2);
Phi_tau_gamma=zeros(Nt,4);
Phi_tau_theta=zeros(Nt,1);

for i=1:Nt
    %price
    S0=S_tau(i,:);
    Ttau=T(i);
    [c1,c1_d1,c1_d2,c1_d1d1,c1_d2d2,c1_d1d2,c1_d2d1,c1_dT]=basket4(S0,Ttau,r,sigma,K1,Rho);
    [c2,c2_d1,c2_d2,c2_d1d1,c2_d2d2,c2_d1d2,c2_d2d1,c2_dT]=basket4(S0,Ttau,r,sigma,K2,Rho);
    [c3,c3_d1,c3_d2,c3_d1d1,c3_d2d2,c3_d1d2,c3_d2d1,c3_dT]=basket4(S0,Ttau,r,sigma,K3,Rho);
    Phi_tau(i)=-2*c2+c1+c3;
    %delta d1
    Phi_tau_delta(i,1)=-2*c2_d1+c1_d1+c3_d1;
    % delta d2
    Phi_tau_delta(i,2)=-2*c2_d2+c1_d2+c3_d2;
    % gamma d11
    Phi_tau_gamma(i,1)=-2*c2_d1d1+c1_d1d1+c3_d1d1;
    Phi_tau_gamma(i,2)=-2*c2_d1d2+c1_d1d2+c3_d1d2;
    Phi_tau_gamma(i,3)=-2*c2_d2d1+c1_d2d1+c3_d2d1;
    Phi_tau_gamma(i,4)=-2*c2_d2d2+c1_d2d2+c3_d2d2;
    % Theta
    Phi_tau_theta(i)=-2*c2_dT + c1_dT+c3_dT;
end




