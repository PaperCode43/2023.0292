function SX=Spath3D(S,t,T,m,J,r,sigma,Rho)
dt=(T-t)/m;
rho=Rho(1,2);
%Z0=normrnd(0,1,m,J);
Z1=normrnd(0,1,m,J);
Z2=rho*Z1+sqrt(1-rho^2)*normrnd(0,1,m,J);
Sm1=S(1)*ones(m+1,J);
Sm2=S(2)*ones(m+1,J);
SX=zeros(m+1,3,J);
for i=1:m
     Sm1(i+1,:)=Sm1(i,:).*exp((r-0.5*sigma(1)^2)*dt+sigma(1)*sqrt(dt)*Z1(i,:));
     Sm2(i+1,:)=Sm2(i,:).*exp((r-0.5*sigma(2)^2)*dt+sigma(2)*sqrt(dt)*Z2(i,:));
end
SX(:,1,:)=Sm1;
SX(:,2,:)=Sm2;
SX(:,3,:)=(t:dt:T)'*ones(1,J);