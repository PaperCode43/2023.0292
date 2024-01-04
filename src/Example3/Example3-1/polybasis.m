function basis=polybasis(X,deg)
basis=ones(size(X,1),1);
for i=1:deg
    basis=[basis,X.^i];
end