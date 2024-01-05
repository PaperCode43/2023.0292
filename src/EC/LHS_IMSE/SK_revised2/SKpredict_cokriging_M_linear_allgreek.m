%%% revised by SHH 2016.2.23 (add MSE calculation)
%%% revised by SHH 2016.3.9 (add hypothesis test, p_Value)

function [ws,w_delta,w_vega,w_rho,w_theta,beta] = SKpredict_cokriging_M_linear_allgreek(model,Xpred,Bpred)
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
R1 = zeros(k,K);R2 = zeros(k,K);R3 = zeros(k,K);R4 = zeros(k,K);
R10 = zeros(k,K);R11 = zeros(k,K);R12 = zeros(k,K);R13 = zeros(k,K);R14 = zeros(k,K);
R20 = zeros(k,K);R21 = zeros(k,K);R22 = zeros(k,K);R23 = zeros(k,K);R24 = zeros(k,K);
R30 = zeros(k,K);R31 = zeros(k,K);R32 = zeros(k,K);R33 = zeros(k,K);R34 = zeros(k,K);
R40 = zeros(k,K);R41 = zeros(k,K);R42 = zeros(k,K);R43 = zeros(k,K);R44 = zeros(k,K);

for i=1:k %Xpred - K x d
    R1(i,:) = 2*theta(1)*(Xpred(:,1)' - X(i,1)) .* R0(i,:);%theta
    R2(i,:) = 2*theta(2)*(Xpred(:,2)' - X(i,2)) .* R0(i,:);%vega
    R3(i,:) = 2*theta(3)*(Xpred(:,3)' - X(i,3)) .* R0(i,:);%rho
    R4(i,:) = 2*theta(4)*(Xpred(:,4)' - X(i,4)) .* R0(i,:);%theta
    R10(i,:) = -R1(i,:);
    R11(i,:) = 2*theta(1)*(1 - 2*theta(1) * (Xpred(:,1)' - X(i,1)).^2).*R0(i,:);
    R12(i,:) = -2*theta(1)*2*theta(2)*(Xpred(:,1)' - X(i,1)).*(Xpred(:,2)' - X(i,2)).* R0(i,:);
    R13(i,:) = -2*theta(1)*2*theta(3)*(Xpred(:,1)' - X(i,1)).*(Xpred(:,3)' - X(i,3)).* R0(i,:);
    R14(i,:) = -2*theta(1)*2*theta(4)*(Xpred(:,1)' - X(i,1)).*(Xpred(:,4)' - X(i,4)).* R0(i,:);
    R20(i,:) = -R2(i,:);
    R21(i,:) = R12(i,:);
    R22(i,:) = 2*theta(2)*(1 - 2*theta(2) * (Xpred(:,2)' - X(i,2)).^2).*R0(i,:);
    R23(i,:) = -2*theta(2)*2*theta(3)*(Xpred(:,2)' - X(i,2)).*(Xpred(:,3)' - X(i,3)).* R0(i,:);
    R24(i,:) = -2*theta(2)*2*theta(4)*(Xpred(:,2)' - X(i,2)).*(Xpred(:,4)' - X(i,4)).* R0(i,:);
    R30(i,:) = -R3(i,:);
    R31(i,:) = R13(i,:);
    R32(i,:) = R23(i,:);
    R33(i,:) = 2*theta(3)*(1 - 2*theta(3) * (Xpred(:,3)' - X(i,3)).^2).*R0(i,:);
    R34(i,:) = -2*theta(3)*2*theta(4)*(Xpred(:,3)' - X(i,3)).*(Xpred(:,4)' - X(i,4)).* R0(i,:);
    R40(i,:) = -R4(i,:);
    R41(i,:) = R14(i,:);
    R42(i,:) = R24(i,:);
    R43(i,:) = R34(i,:);
    R44(i,:) = 2*theta(4)*(1 - 2*theta(4) * (Xpred(:,4)' - X(i,4)).^2).*R0(i,:);
end


Rspred=[R0;R1;R2;R3;R4];
Rdeltapred=[R10;R11;R12;R13;R14];
Rvegapred=[R20;R21;R22;R23;R24];
Rrhopred=[R30;R31;R32;R33;R34];
Rthetapred=[R40;R41;R42;R43;R44];

% calculate responses at prediction points 
%f = Bpred*beta + Rpred'/Sigma*Z;
e1=ones(5*k,1);
% ws=((Rspred + e1 * (1-e1'*((Sigmaplus)\Rspred))/(e1'*((Sigmaplus)\e1)))'/(Sigmaplus))';
% w_delta=((Rdeltapred + e1 * (0-e1'*((Sigmaplus)\Rdeltapred))/(e1'*((Sigmaplus)\e1)))'/(Sigmaplus))';%% previous using 1, which is wrong
% w_vega=((Rvegapred + e1 * (0-e1'*((Sigmaplus)\Rvegapred))/(e1'*((Sigmaplus)\e1)))'/(Sigmaplus))';
% w_rho=((Rrhopred + e1 * (0-e1'*((Sigmaplus)\Rrhopred))/(e1'*((Sigmaplus)\e1)))'/(Sigmaplus))';
% w_theta=((Rthetapred + e1 * (0-e1'*((Sigmaplus)\Rthetapred))/(e1'*((Sigmaplus)\e1)))'/(Sigmaplus))';
ws=((Rspred)'/(Sigmaplus))';
w_delta=((Rdeltapred)'/(Sigmaplus))';%% previous using 1, which is wrong
w_vega=((Rvegapred)'/(Sigmaplus))';
w_rho=((Rrhopred )'/(Sigmaplus))';
w_theta=((Rthetapred )'/(Sigmaplus))';

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




