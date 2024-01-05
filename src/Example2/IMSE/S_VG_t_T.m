function [S_VG_t,S_VG_T]=S_VG_t_T(S_0,t,T,n,N,r_a,r_n,sigma_vg,nu_vg,theta_vg)
%n: # time discrete number
%N: # sample size
phi0=-1/nu_vg*log(1-theta_vg*nu_vg-0.5*sigma_vg^2*nu_vg);%martingle component
dt=T/n;
if mod(t,dt)~=0
    error('t/dt is not integer')
end
nt=round(t/dt);
S_VG_t=S_0*ones(N,1);
for i=1:nt
    VG_p=VGp(sigma_vg,nu_vg,theta_vg,N,dt);
    S_VG_t=S_VG_t.*exp((r_a-phi0)*dt + VG_p');
    %S_VG_t=S_VG_t+( VG_p');
end
S_VG_T=S_VG_t;
for i=nt+1:n
    VG_p=VGp(sigma_vg,nu_vg,theta_vg,N,dt);
    S_VG_T=S_VG_T.*exp((r_n-phi0)*dt + VG_p');
end
    
    