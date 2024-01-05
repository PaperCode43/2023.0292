%% this file to calculate the analytical formula of portfolio value
function [Phi_tau,Phi_tau_delta,Phi_tau_sigma,Phi_tau_theta]=BS_price(S_tau,r,sigma,T,K1)
Nt=size(S_tau,1);
Phi_tau=zeros(Nt,1);
Phi_tau_delta=zeros(Nt,1);
Phi_tau_sigma=zeros(Nt,1);
Phi_tau_theta=zeros(Nt,1);

for i=1:Nt
    %price
    S0=S_tau(i);
    Ttau=T(i);
    [c1,c1_d1,c1_dsigma,c1_dT]=bls_Euro_call(S0,Ttau,r,sigma,K1);
    Phi_tau(i)=c1;
    %delta
    Phi_tau_delta(i)=c1_d1;
    % sigma
    Phi_tau_sigma(i)=c1_dsigma;
    % Theta
    Phi_tau_theta(i)=c1_dT;
end




