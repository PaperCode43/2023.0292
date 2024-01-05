function SX=Spath2D(S,t,T,m,J,r,sigma)
dt=(T-t)/m;
Z1=normrnd(0,1,m,J);
Sm=S*ones(m+1,J);
SX=zeros(m+1,2,J);
for i=1:m
     Sm(i+1,:)=Sm(i,:).*exp((r-0.5*sigma^2)*dt+sigma*sqrt(dt)*Z1(i,:));
end
SX(:,1,:)=Sm;
SX(:,2,:)=(t:dt:T)'*ones(1,J);