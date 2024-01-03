%%% revised by SHH 2016.2.23 (add MSE calculation)
%%% revised by SHH 2016.3.9 (add hypothesis test, p_Value)

function [ws,wd,beta] = SKpredict_cokriging_M_linear(model,Xpred,Bpred)
% make predictions at prediction points using a stochastic kriging model  
% model = output of SKfit
% Xpred = (K x d) matrix of prediction points
% Bpred = (K x b) matrix of basis functions at each prediction point
%         The first column must be a column of ones!
% f = (K x 1) predictions at predictions points
% 
% Exmaples
%      SK_gau  = SKpredict(skriging_model,XK,ones(K,1));
% Based on parameter estimates skriging_model obtained from SKfit.m,
% use SK model to predict the values at prediction points XK with constant
% prediction trend, Bpred = ones(K,1)

% retrieve model parameters from model structure obtained from SKfit
X = model.X;
% minX = model.minX;
% maxX = model.maxX;
[k, d] = size(X);
theta = model.theta;
%theta=0.1;
gammaP = model.gamma;
beta = model.beta;
%beta=0;
Sigmaplus = model.Sigmaplus;
Z = model.Z;
tau2 = model.tausquared;
%tau2=0.5;

% simple check for dimensions of Xpred and X
K = size(Xpred,1);     % number of prediction points
if (size(Xpred,2)~=d)
    error('Prediction points and design points must have the same dimension (number of columns).');
end
if (size(Bpred,1)~=K)
    error('Basis function and prediction point matrices must have the same number of rows.')
end
if not(all(Bpred(:,1)==1))
    error('The first column of the basis function matrix must be ones.')
end

% calculate distance matrix for prediction points
% Xpred = (Xpred - repmat(minX,K,1)) ./ repmat(maxX-minX,K,1);
if gammaP == 2
    distXpred =  abs(repmat(reshape(Xpred', [1 d K]),[k,1,1]) ...
        - repmat(X,[1 1 K])).^2;
else
    distXpred =  abs(repmat(reshape(Xpred', [1 d K]),[k,1,1]) ...
        - repmat(X,[1 1 K]));
end

% calculate correlations between prediction points and design points
D = distXpred;
if gammaP == 3
    T = repmat(reshape(theta,[1 d 1]),[k 1 K]);
    Rspred = tau2*prod(((D<=(T./2)).*(1-6*(D./T).^2+6*(D./T).^3) ...
        +((T./2)<D & D<=T).*(2*(1-D./T).^3)),2);
else
    Rspred = tau2*exp(sum(-D.*repmat(reshape(theta,[1 d]),[k 1 K]),2));
end
R0 = reshape(Rspred,[k K 1]); %k x K
R1 = zeros(k,K);
R10 = zeros(k,K);
R11 = zeros(k,K);
for i=1:k %Xpred - K x d
    R1(i,:) = 2*theta(1)*(Xpred(:,1)' - X(i,1)) .* R0(i,:);
    R10(i,:) = -R1(i,:);
    R11(i,:) = 2*theta(1)*(1 - 2*theta(1) * (Xpred(:,1)' - X(i,1)).^2).*R0(i,:);
end


Rspred=[R0;R1];
Rdpred=[R10;R11];

% calculate responses at prediction points 
%f = Bpred*beta + Rpred'/Sigma*Z;
e1=ones(2*k,1);
%ws=((Rspred + e1 * (1-e1'*((Sigmaplus)\Rspred))/(e1'*((Sigmaplus)\e1)))'/(Sigmaplus))';
%wd=((Rdpred + e1 * (0-e1'*((Sigmaplus)\Rdpred))/(e1'*((Sigmaplus)\e1)))'/(Sigmaplus))';
ws=((Rspred)'/(Sigmaplus))';
wd=((Rdpred)'/(Sigmaplus))';
%fun=@(x)fobjc(x,Sigma,Rpred);
x0=ones(k,1)/k;
%fobjc(x0,Sigma,Rpred)
%w = fmincon(fun,w0,e1',1,[],[],zeros(k,1),[]);


function [f,g] = fobjc(w,Sigma,Rpred)
f = -(-0.5*w'*Sigma*w+Rpred'*w);
g = -(-Sigma*w+Rpred);



%%%%%%
% if nargout > 1        %%% calculate the MSE _ add by SHH
%     MSE = zeros(size(Xpred),1);
%     B = model.B;
%     i = 1;
%     for X0 = Xpred'
%         r = Rpred(:,i);      
%         delta = Bpred(i,:)' - B'/Sigma*r;    
%         MSE(i,1) = tau2 - r'/Sigma*r + delta'*(B'/Sigma*B)^(-1)*delta;      
%         i=i+1;
%     end
% end




