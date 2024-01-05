function [S_VG,G,Z]=S_VGT(S_0,T,n,N,r_n,sigma_vg,nu_vg,theta_vg)
%n: # time discrete number
%N: # sample size
phi0=-1/nu_vg*log(1-theta_vg*nu_vg-0.5*sigma_vg^2*nu_vg);%martingle component
dt=T/n;
S_VG_t=S_0*ones(N,1);
VG_p=zeros(N,n);
G=zeros(N,n);
Z=zeros(N,n);
S_VG=zeros(N,n);
for i=1:n
    [VG_p(:,i),G(:,i),Z(:,i)]=VGp(sigma_vg,nu_vg,theta_vg,N,dt);
    S_VG_t=S_VG_t.*exp((r_n-phi0)*dt + VG_p(:,i));
    S_VG(:,i) = S_VG_t;
end

function [VG,G,Z]=VGp(sigma,nu,theta,N,T)
G=gamrnd(T/nu,nu,N,1);
Z=randn(N,1);
WG=Z.*sqrt(G);
VG=theta*G+sigma*WG;

    
    