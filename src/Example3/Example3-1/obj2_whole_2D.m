function objstage2_integral=obj2_whole_2D(par,X,Xtest,Vhat,D0,gammaP,Y,beta,tau2,r,sigma,T,m,J,Vmatrix)

[n,~]=size(Xtest);
obj_temp=zeros(n,m);
%rng('default');
for i=1:n
    SX=Spath2D(Xtest(i,1),Xtest(i,2),T,m,J,r,sigma);
    for j=1:m
            A=reshape(SX(j,:,:),2,J);
            obj_temp(i,j)=exp(-r*(T-Xtest(i,2)))*(obj2_nonoise(par,X,A',Vhat,D0,gammaP,Y,beta,tau2,r,sigma,Vmatrix))*(T-Xtest(i,2)/m);
    end
end
objstage2_integral=mean((sum(obj_temp,2)).^2);