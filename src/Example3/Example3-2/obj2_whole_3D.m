function objstage2_integral=obj2_whole_3D(par,X,Xtest,Vhat,D0,gammaP,Y,beta,tau2,r,sigma,Rho,T,m,J,Vmatrix)

[n,~]=size(Xtest);
obj_temp=zeros(n,m);
%rng('default');
for i=1:n
    SX=Spath3D(Xtest(i,1:2),Xtest(i,3),T,m,J,r,sigma,Rho);
    for j=1:m
            A=reshape(SX(j,:,:),3,J);
            obj_temp(i,j)=exp(-r*(T-Xtest(i,3)))*(obj2_nonoise2(par,X,A',Vhat,D0,gammaP,Y,beta,tau2,r,sigma,Rho,Vmatrix))*(T-Xtest(i,3)/m);
    end
end
objstage2_integral=mean((sum(obj_temp,2)).^2);