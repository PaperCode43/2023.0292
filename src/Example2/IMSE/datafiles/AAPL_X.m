function X=AAPL_X(S0,sigma,theta)
u=lhsdesign(30,3);
%%%%X=(S0,sigma,r,T)
X=zeros(30,3);
X(:,1) = (1+(2*u(:,1)-1)*0.3)*S0;
X(:,2) = (1+(2*u(:,2)-1)*0.3)*sigma;
X(:,3) = (1+(2*u(:,3)-1)*0.3)*theta;