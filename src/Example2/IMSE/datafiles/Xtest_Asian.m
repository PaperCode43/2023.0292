%%Apple 
%%% [S0,sigma,theta]
function [Xtest,Y_true_test]=Xtest_Asian(S0,sigma_vg,theta_vg,K,T,n,r_n,nu_vg,N)
p=sobolset(3);
u=p(2:101,:);
Xtest=zeros(100,3);
Xtest(:,1) = (1+(2*u(:,1)-1)*0.31)*S0;
Xtest(:,2) = (1+(2*u(:,2)-1)*0.31)*sigma_vg;
Xtest(:,3) = (1+(2*u(:,3)-1)*0.31)*theta_vg;
Y_true_test.V=zeros(100,1);
Y_true_test.dS0=zeros(100,1);
Y_true_test.dsigma=zeros(100,1);
Y_true_test.dtheta=zeros(100,1);
for i=1:100
    
    [Y_true_test.V(i),Y_true_test.dS0(i),Y_true_test.dsigma(i),Y_true_test.dtheta(i),~] = VG_call_Asian(Xtest(i,1),K,T,n,N,r_n,Xtest(i,2),nu_vg,Xtest(i,3));
end
