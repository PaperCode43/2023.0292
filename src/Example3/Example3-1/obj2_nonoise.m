function objstage2=obj2_nonoise(par,X,Xtest,Vhat,D0,gammaP,Y,beta,tau2,r,sigma,Vmatrix)
theta=par;
[k, d] = size(X);
K = size(Xtest,1);
Rhat = corrfun(theta,D0,gammaP);
B=ones(k,1);
Sigmahat_s = tau2*Rhat + diag(Vhat);
temp1 = B'/Sigmahat_s;
temp2 = temp1*B;
temp3 = temp1*Y(1:k);
beta = temp2\temp3;


C00 = tau2*Rhat;
C01=zeros(k,k);
C11=zeros(k,k);
%NOTICE here we specitic to 2-d
for i=1:k
    for h=1:k
        C01(i,h) = 2*theta(1)*(X(i,1)-X(h,1))*C00(i,h);
        C11(i,h) = 2*theta(1)*(1 - 2*theta(1)*(X(i,1)-X(h,1))^2)*C00(i,h);
    end
end
C10=-C01;
Sigma_M=[C00,C01;C10,C11];
Sigmahat = Sigma_M + Vmatrix;


if gammaP == 2
    distXpred =  abs(repmat(reshape(Xtest', [1 d K]),[k,1,1]) ...
        - repmat(X,[1 1 K])).^2;
else
    distXpred =  abs(repmat(reshape(Xtest', [1 d K]),[k,1,1]) ...
        - repmat(X,[1 1 K]));
end

% calculate correlations between prediction points and design points
D = distXpred;
if gammaP == 3
    T = repmat(reshape(theta,[1 d 1]),[k 1 K]);
    Rspred = tau2*prod(((D<=(T./2)).*(1-6*(D./T).^2+6*(D./T).^3) ...
        +((T./2)<D & D<=T).*(2*(1-D./T).^3)),2);
else
    Rspred = tau2*exp(sum(-D.*repmat(reshape(theta,[1 d]),[k 1 K]),2));
end
R0 = reshape(Rspred,[k K 1]); %k x K
R1 = -2*theta(1)*(ones(k,1)*(Xtest(:,1)') - X(:,1)*ones(1,K)).* R0;
R2 = -2*theta(2)*(ones(k,1)*(Xtest(:,2)') - X(:,2)*ones(1,K)).* R0;
R11 = 2*theta(1)*(2*theta(1)*(ones(k,1)*(Xtest(:,1)') - X(:,1)*ones(1,K)).^2-1).* R0;
%R22 = 2*theta(2)*(2*theta(2)*(ones(k,1)*(Xtest(:,2)') - X(:,2)*ones(1,K)).^2-1).* R0;
R12 = 4*theta(1)*theta(2)*((ones(k,1)*(Xtest(:,1)') - X(:,1)*ones(1,K)).*(ones(k,1)*(Xtest(:,2)') - X(:,2)*ones(1,K))).* R0;
R111 = (12*theta(1)^2*(ones(k,1)*(Xtest(:,1)') - X(:,1)*ones(1,K)) - 8*theta(1)^3*(ones(k,1)*(Xtest(:,1)') - X(:,1)*ones(1,K)).^3).* R0;
%R112 = (-2*theta(2)*(ones(k,1)*(Xtest(:,2)') - X(:,2)*ones(1,K))).*R11;
betam=[ones(k,1)*beta;zeros(k,1)];

tem=inv(Sigmahat)*(Y - betam);
ob1 = [R2;R12]'*tem;
ob2 = r*Xtest(:,1).*([R1;R11]'*tem);
ob3 = 0.5*sigma^2*Xtest(:,1).^2.*([R11;R111]'*tem);
ob4=(beta+[R0;R1]'*tem);
ob=(ob1+ob2+ob3-r*ob4);
objstage2=sum(ob);

% mob2=mean(ob2.^2)
% mob3=mean(ob3.^2)
% mob4=mean(ob4.^2)
% cdn=cond(Sigmahat)