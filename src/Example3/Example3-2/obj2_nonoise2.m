function objstage2=obj2_nonoise2(par,X,Xtest,Vhat,D0,gammaP,Y,beta,tau2,r,sigma,Rho,Vmatrix)
theta=par;
rho=Rho(1,2);
[k, d] = size(X);
K = size(Xtest,1);
Rhat = corrfun(theta,D0,gammaP);

%%revised by Jiang 20230208
B=ones(k,1);
Sigmahat_s = tau2*Rhat + diag(Vhat);
temp1 = B'/Sigmahat_s;
temp2 = temp1*B;
temp3 = temp1*Y(1:k);
beta = temp2\temp3;


C00 = tau2*Rhat;
C01=zeros(k,k);C02=zeros(k,k);
C11=zeros(k,k);C12=zeros(k,k);
C22=zeros(k,k);
%NOTICE here we specitic to 2-d
for i=1:k
    for h=1:k
        C01(i,h) = 2*theta(1)*(X(i,1)-X(h,1))*C00(i,h);
        C02(i,h) = 2*theta(2)*(X(i,2)-X(h,2))*C00(i,h);
        C11(i,h) = 2*theta(1)*(1 - 2*theta(1)*(X(i,1)-X(h,1))^2)*C00(i,h);
        C12(i,h) = -2*theta(1)*(X(i,1)-X(h,1))*2*theta(2)*(X(i,2)-X(h,2))*C00(i,h);
        C22(i,h) = 2*theta(2)*(1 - 2*theta(2)*(X(i,2)-X(h,2))^2)*C00(i,h);
    end
end
C10=-C01;C20=-C02;C21=C12;
Sigma_M=[C00,C01,C02;C10,C11,C12;C20,C21,C22];
Sigmahat = Sigma_M + 0.00001*Vmatrix;



%Sigmahat  = tau2*Rhat;

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
R3 = -2*theta(3)*(ones(k,1)*(Xtest(:,3)') - X(:,3)*ones(1,K)).* R0;
R11 = 2*theta(1)*(2*theta(1)*(ones(k,1)*(Xtest(:,1)') - X(:,1)*ones(1,K)).^2-1).* R0;
R22 = 2*theta(2)*(2*theta(2)*(ones(k,1)*(Xtest(:,2)') - X(:,2)*ones(1,K)).^2-1).* R0;
R12 = 4*theta(1)*theta(2)*((ones(k,1)*(Xtest(:,1)') - X(:,1)*ones(1,K)).*(ones(k,1)*(Xtest(:,2)') - X(:,2)*ones(1,K))).* R0;
R13 = 4*theta(1)*theta(3)*((ones(k,1)*(Xtest(:,1)') - X(:,1)*ones(1,K)).*(ones(k,1)*(Xtest(:,3)') - X(:,3)*ones(1,K))).* R0;
R23 = 4*theta(2)*theta(3)*((ones(k,1)*(Xtest(:,2)') - X(:,2)*ones(1,K)).*(ones(k,1)*(Xtest(:,3)') - X(:,3)*ones(1,K))).* R0;
R111 = (12*theta(1)^2*(ones(k,1)*(Xtest(:,1)') - X(:,1)*ones(1,K)) - 8*theta(1)^3*(ones(k,1)*(Xtest(:,1)') - X(:,1)*ones(1,K)).^3).* R0;
R222 = (12*theta(2)^2*(ones(k,1)*(Xtest(:,2)') - X(:,2)*ones(1,K)) - 8*theta(2)^3*(ones(k,1)*(Xtest(:,2)') - X(:,2)*ones(1,K)).^3).* R0;
R112 = (-2*theta(2)*(ones(k,1)*(Xtest(:,2)') - X(:,2)*ones(1,K))).*R11;
R221 = (-2*theta(1)*(ones(k,1)*(Xtest(:,1)') - X(:,1)*ones(1,K))).*R22;

betam=[ones(k,1)*beta;zeros(2*k,1)];
%kap1=[R3;R13;R23]';
tem=inv(Sigmahat)*(Y - betam);

%ob1 = [R3;R13;R23]'/(Sigmahat)*(Y - betam);
ob1 = [R3;R13;R23]'*tem;
%kap2=(r*(ones(3*k,1)*Xtest(:,1)').*[R1;R11;R12])'+(r*(ones(3*k,1)*Xtest(:,2)').*[R2;R12;R22])';
%ob2=r*Xtest(:,1).*([R1;R11;R12]'/(Sigmahat)*(Y - betam))+r*Xtest(:,2).*([R2;R12;R22]'/(Sigmahat)*(Y - betam));
ob2=r*Xtest(:,1).*([R1;R11;R12]'*tem)+r*Xtest(:,2).*([R2;R12;R22]'*tem);
%kap3=((0.5*sigma(1)^2*(ones(3*k,1)*Xtest(:,1)').^2.*[R11;R111;R112])+(0.5*sigma(2)^2*(ones(3*k,1)*Xtest(:,2)').^2.*[R22;R221;R222])+...
 %   (sigma(1)*sigma(2)*rho*(ones(3*k,1)*Xtest(:,1)').*(ones(3*k,1)*Xtest(:,2)').*[R12;R112;R221]))';
%ob3 = 0.5*sigma(1)^2*Xtest(:,1).^2.*([R11;R111;R112]'/(Sigmahat)*(Y - betam)) + 0.5*sigma(2)^2*Xtest(:,2).^2.*([R22;R221;R222]'/(Sigmahat)*(Y - betam))...
   % + sigma(1)*sigma(2)*rho*Xtest(:,1).*Xtest(:,2).*([R12;R112;R221]'/(Sigmahat)*(Y - betam));
ob3 = 0.5*sigma(1)^2*Xtest(:,1).^2.*([R11;R111;R112]'*tem) + 0.5*sigma(2)^2*Xtest(:,2).^2.*([R22;R221;R222]'*tem)...
    + sigma(1)*sigma(2)*rho*Xtest(:,1).*Xtest(:,2).*([R12;R112;R221]'*tem);
%kap4=[R0;R1;R2]';

%ob4=(beta+[R0;R1;R2]'/(Sigmahat)*(Y - betam));
ob4=(beta+[R0;R1;R2]'*tem);

%ob1=kap1/(Sigmahat)*(Y - betam);
%ob2=(kap2/(Sigmahat)*(Y - betam));
%ob3=(kap3/(Sigmahat)*(Y - betam));
%ob4=(beta+kap4/(Sigmahat)*(Y - betam));
ob=(ob1+ob2+ob3-r*ob4);

objstage2=sum(ob);
