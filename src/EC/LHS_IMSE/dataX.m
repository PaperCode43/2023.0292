p=sobolset(4);
u=p(2:21,:);
%%%%X=(S0,sigma,r,T)
X=zeros(20,4);
X(:,1) = (1+(2*u(:,1)-1)*0.2)*100;
X(:,2) = u(:,2)*0.3;
X(:,3) = u(:,3)*0.1;
X(:,4) = u(:,4)*2;