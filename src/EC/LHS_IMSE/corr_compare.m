%% Correlation comparison
function [c1,c2,c_mM,c_Mm,d1,d2,d_mM,d_Mm]=corr_compare(n,d,rep)

X_mini = readmatrix('miniLHD.xlsx');
X_Maxi = readmatrix('MmLHD.xlsx');

c1=zeros(rep,1);c2=zeros(rep,1);c_mM=zeros(rep,1);c_Mm=zeros(rep,1);
d1=zeros(rep,1);d2=zeros(rep,1);d_mM=zeros(rep,1);d_Mm=zeros(rep,1);
for i=1:rep
    rng(i)
    X1 = lhsdesign(n,d);
    y1 = corr(X1);
    c1(i) = (sum(y1(:).^2) - d)/2;
    d1(i)=mindistance(X1);

    rng(i)
    X2 = lhsdesign(10,d,'Criterion','correlation');
    y2 = corr(X2);
    c2(i) = (sum(y2(:).^2) - d)/2;
    d2(i)=mindistance(X2);

    X_mM = X_mini(1+(i-1)*n:i*n,1:d);
    y_mM = corr(X_mM);
    c_mM(i) = (sum(y_mM(:).^2) - d)/2;
    d_mM(i)=mindistance(X_mM);

    X_Mm = X_Maxi(1+(i-1)*n:i*n,1:d);
    y_Mm = corr(X_Mm);
    c_Mm(i) = (sum(y_Mm(:).^2) - d)/2;
    d_Mm(i)=mindistance(X_Mm);
end

% figure(1)
% hold on
% boxplot(([c1,c2,c_mM,c_Mm]),...
%     'Labels',{'LHD-Matlab','LHD-Matlab-corr','LHD-mM','LHD-Mm'})
% ylabel('Averaged Correlation')
% set(gca,'FontSize',16);
% set(findobj(gca,'Type','text'),'FontSize',16)
% hold off
% 
% figure(2)
% hold on
% boxplot(([d1,d2,d_mM,d_Mm]),...
%     'Labels',{'LHD-Matlab','LHD-Matlab-corr','LHD-mM','LHD-Mm'})
% ylabel('Minimum Distance')
% set(gca,'FontSize',16);
% set(findobj(gca,'Type','text'),'FontSize',16)
% hold off
md1=mean(c1)
std1=std(c1)
md2=mean(c2)
std2=std(c2)

function d=mindistance(X)
n=size(X,1);
id=1;
dm=zeros(n*(n-1)/2,1);
for i=1:n
    for j=1:i-1
        dm(id)=sqrt(sum((X(i,:) - X(j,:)).^2));
        id=id+1;
    end
end
d=min(dm);
        