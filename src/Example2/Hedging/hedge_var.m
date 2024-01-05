%% for the left panel, please run "rng(1);[p1,p2,S]=Portfolio_gap(63);"
rep=40;
Nday=63;
portfolio_ESK=zeros(rep,Nday-1);
portfolio_SK=zeros(rep,Nday-1);
SS=zeros(Nday+1,5,rep);
parfor i=1:rep
    rng(i);
[p1,p2,S]=Portfolio_gap(Nday);
portfolio_ESK(i,:)=p1;
portfolio_SK(i,:)=p2;
SS(:,:,i)=S;
end

var_port_ESK = var(portfolio_ESK);
var_port_SK = var(portfolio_SK);
plot(1:Nday-1,var_port_ESK,'r*-',1:Nday-1,var_port_SK,'bo--','LineWidth',2)
save('63days_var')
