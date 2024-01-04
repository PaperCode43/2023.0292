function [w_delta,w_gamma,w_theta,beta] = SKpredict_M_linear_der_delta_gamma_theta_pde(model,Xpred,Bpred)
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
Sigma = model.Sigma;
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
    Rpred = tau2*prod(((D<=(T./2)).*(1-6*(D./T).^2+6*(D./T).^3) ...
        +((T./2)<D & D<=T).*(2*(1-D./T).^3)),2);
else
    Rpred = tau2*exp(sum(-D.*repmat(reshape(theta,[1 d]),[k 1 K]),2));
end
R0 = reshape(Rpred,[k K 1]);
R_delta = zeros(k,K);
R_gamma = zeros(k,K);
R_theta = zeros(k,K);

for i=1:k %Xpred - K x d
     R_delta(i,:) = -2*theta(1)*(Xpred(:,1)' - X(i,1)).* R0(i,:);
     R_theta(i,:)  = -2*theta(2)*(Xpred(:,2)' - X(i,2)).* R0(i,:);
     R_gamma(i,:) = 2*theta(1)*(2*theta(1)*(Xpred(:,1)' - X(i,1)).^2-1).* R0(i,:);
end

% calculate responses at prediction points 
w_delta=(R_delta'/(Sigma))';
w_theta=(R_theta'/(Sigma))';
w_gamma=(R_gamma'/(Sigma))';







