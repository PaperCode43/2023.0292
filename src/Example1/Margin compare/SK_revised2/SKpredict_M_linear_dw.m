%%% revised by SHH 2016.2.23 (add MSE calculation)
%%% revised by SHH 2016.3.9 (add hypothesis test, p_Value)

function dw = SKpredict_M_linear_dw(model,Xpred,Bpred,idx)
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
if idx>d
    error('exceeding the dimension of x')
end
theta = model.theta;
%theta=0.1;
gammaP = model.gamma;
%beta=0;
Sigma = model.Sigma;
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
    Rpred = tau2*prod(((D<=(T./2)).*(1-6*(D./T).^2+6*(D./T).^3) ...
        +((T./2)<D & D<=T).*(2*(1-D./T).^3)),2);
else
    Rpred = tau2*exp(sum(-D.*repmat(reshape(theta,[1 d]),[k 1 K]),2));
end
Rpred = reshape(Rpred,[k K 1]);
theta_i=theta(idx);
if gammaP == 3
    error('We DO NOT consider case: gamma =3')
elseif gammaP == 2
    term1 = -2*theta_i*(ones(k,1)*(Xpred(:,idx))'-X(:,idx)*ones(1,K));
elseif gammaP == 1
    term1 = -theta_i.*(ones(k,1)*(Xpred(:,idx))'>X(:,idx)*ones(1,K));
end
dRpred = term1.*Rpred;

% calculate responses at prediction points 
%f = Bpred*beta + Rpred'/Sigma*Z;
e1=ones(k,1);
dw=((dRpred + e1 * (0-e1'*inv(Sigma)*dRpred)/(e1'*inv(Sigma)*e1))'*inv(Sigma))';

%fobjc(x0,Sigma,Rpred)
%w = fmincon(fun,w0,e1',1,[],[],zeros(k,1),[]);






