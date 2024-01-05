%%% revised by SHH 2016.2.23 (add MSE calculation)
%%% revised by SHH 2016.3.9 (add hypothesis test, p_Value)

function [ws,w_dS01,w_dS02,w_gamma11,w_gamma12,w_gamma22,w_dt,beta] = SKpredict_cokriging_M_linear_delta_pde2(model,Xpred,Bpred,Vmatrix)
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
%Sigmaplus = model.Sigma;
Z = model.Z;
tau2 = model.tausquared;

%%
% % Normalize data by scaling each dimension from 0 to 1
% minX = min(X);  
% maxX = max(X);
% X = (X - repmat(minX,k,1)) ./ repmat(maxX-minX,k,1);

% calculate the distance matrix between points (copied from DACE)
% distances are recorded for each dimension separately
ndistX = k*(k-1) / 2;        % max number of non-zero distances
ijdistX = zeros(ndistX, 2);  % initialize matrix with indices
distX = zeros(ndistX, d);    % initialize matrix with distances
temp = 0;
for i = 1 : k-1
    temp = temp(end) + (1 : k-i);
    ijdistX(temp,:) = [repmat(i, k-i, 1) (i+1 : k)']; 
    distX(temp,:) = repmat(X(i,:), k-i, 1) - X(i+1:k,:); 
end
IdistX = sub2ind([k k],ijdistX(:,1),ijdistX(:,2));

% distance matrix raised to the power of gamma
D0 = zeros(k,k,d);
for p=1:d
    temp = zeros(k);
    if (gammaP == 3)
        temp(IdistX) = abs(distX(:,p));
    else
        temp(IdistX) = -abs(distX(:,p)).^gammaP;
    end
    D0(:,:,p) = temp+temp';
end

Rhat = corrfun(theta,D0,gammaP);
C00 = tau2*Rhat;
C01=zeros(k,k);C02=zeros(k,k);
C11=zeros(k,k);C12=zeros(k,k);
C22=zeros(k,k);
%NOTICE here we specitic to 2-d
for i=1:k
    for h=1:k
        C01(i,h) = 2*theta(1)*(X(i,1)-X(h,1))*C00(i,h);
        C02(i,h) = 2*theta(2)*(X(i,2)-X(h,2))*C00(i,h);
        C11(i,h) = 2*theta(1)*(1 - 2*theta(1)*(X(i,1)-X(h,1))^2)*C00(i,h);
        C12(i,h) = -2*theta(1)*(X(i,1)-X(h,1))*2*theta(2)*(X(i,2)-X(h,2))*C00(i,h);
        C22(i,h) = 2*theta(2)*(1 - 2*theta(2)*(X(i,2)-X(h,2))^2)*C00(i,h);
    end
end
C10=-C01;C20=-C02;C21=C12;
Sigma_M=[C00,C01,C02;C10,C11,C12;C20,C21,C22];
Sigmaplus = Sigma_M + Vmatrix;


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
R1 = zeros(k,K);R2 = zeros(k,K);
R10 = zeros(k,K);R11 = zeros(k,K);R12 = zeros(k,K);
R20 = zeros(k,K);R21 = zeros(k,K);R22 = zeros(k,K);
Rt1 =zeros(k,K);Rt2 =zeros(k,K);Rt3 =zeros(k,K);
Rgamma11_1 =zeros(k,K); Rgamma11_2 =zeros(k,K); Rgamma11_3 =zeros(k,K); 
Rgamma22_1 =zeros(k,K); Rgamma22_2 =zeros(k,K); Rgamma22_3 =zeros(k,K); 
Rgamma12_1 =zeros(k,K); Rgamma12_2 =zeros(k,K); Rgamma12_3 =zeros(k,K);



for i=1:k %Xpred - K x d
    R1(i,:) = 2*theta(1)*(Xpred(:,1)' - X(i,1)) .* R0(i,:);%dS01
    R2(i,:) = 2*theta(2)*(Xpred(:,2)' - X(i,2)) .* R0(i,:);%dS02
    
    R10(i,:) = -R1(i,:);
    R11(i,:) = 2*theta(1)*(1 - 2*theta(1) * (Xpred(:,1)' - X(i,1)).^2).*R0(i,:);
    R12(i,:) = -2*theta(1)*2*theta(2)*(Xpred(:,1)' - X(i,1)).*(Xpred(:,2)' - X(i,2)).* R0(i,:);
    R20(i,:) = -R2(i,:);
    R21(i,:) = R12(i,:);
    R22(i,:) = 2*theta(2)*(1 - 2*theta(2) * (Xpred(:,2)' - X(i,2)).^2).*R0(i,:);
    
    %taking derivative on w to obtain other weights of derivatives (theta on surface, gamma11 on delta1, gamma22 on delta2, gamma12 on delta1)
    Rt1(i,:)=-2*theta(3)*(Xpred(:,3)' - X(i,3)).* R0(i,:);
    Rt2(i,:)= - 4*theta(1)*theta(3)*(Xpred(:,1)' - X(i,1)).*(Xpred(:,3)' - X(i,3)).* R0(i,:);
    Rt3(i,:)= - 4*theta(2)*theta(3)*(Xpred(:,2)' - X(i,2)).*(Xpred(:,3)' - X(i,3)).* R0(i,:);
    Rgamma11_1(i,:) = -R11(i,:);
    Rgamma11_2(i,:) = (-4*theta(1)^2*(Xpred(:,1)' - X(i,1)) + 8*theta(1)^3*(Xpred(:,1)' - X(i,1)).^3 - 8*theta(1)^2*(Xpred(:,1)' - X(i,1))).*R0(i,:);
    Rgamma11_3(i,:) = (-4*theta(1)*theta(2)*(Xpred(:,2)' - X(i,2)) + 8*theta(1)^2*theta(2)*(Xpred(:,1)' - X(i,1)).^2.*(Xpred(:,2)' - X(i,2))).*R0(i,:);   
    Rgamma22_1(i,:) = -R22(i,:);
    Rgamma22_2(i,:) = (-4*theta(1)*theta(2)*(Xpred(:,1)' - X(i,1)) + 8*theta(1)*theta(2)^2*(Xpred(:,1)' - X(i,1)).*(Xpred(:,2)' - X(i,2)).^2).*R0(i,:);
    Rgamma22_3(i,:) = (-4*theta(2)^2*(Xpred(:,2)' - X(i,2)) + 8*theta(2)^3*(Xpred(:,2)' - X(i,2)).^3 - 8*theta(2)^2*(Xpred(:,2)' - X(i,2))).*R0(i,:);
    Rgamma12_1(i,:) = -R12(i,:);
    Rgamma12_2(i,:) = Rgamma11_3(i,:);
    Rgamma12_3(i,:) = Rgamma22_2(i,:);
    
    
end


Rspred=[R0;R1;R2];
RdS01pred=[R10;R11;R12];
RdS02pred=[R20;R21;R22];
Rdtpred=[Rt1;Rt2;Rt3];
Rgamma11=[Rgamma11_1;Rgamma11_2;Rgamma11_3];
Rgamma22=[Rgamma22_1;Rgamma22_2;Rgamma22_3];
Rgamma12=[Rgamma12_1;Rgamma12_2;Rgamma12_3];


% calculate responses at prediction points 
%f = Bpred*beta + Rpred'/Sigma*Z;
% ws=((Rspred + e1 * (1-e1'*((Sigmaplus)\Rspred))/(e1'*((Sigmaplus)\e1)))'/(Sigmaplus))';
% w_dS0=((RdS0pred + e1 * (0-e1'*((Sigmaplus)\RdS0pred))/(e1'*((Sigmaplus)\e1)))'/(Sigmaplus))';%% previous using 1, which is wrong
% w_dsigma=((Rdsigmapred + e1 * (0-e1'*((Sigmaplus)\Rdsigmapred))/(e1'*((Sigmaplus)\e1)))'/(Sigmaplus))';
% w_dtheta=((Rdthetapred + e1 * (0-e1'*((Sigmaplus)\Rdthetapred))/(e1'*((Sigmaplus)\e1)))'/(Sigmaplus))';
ws=((Rspred)'/(Sigmaplus))';
w_dS01=((RdS01pred)'/(Sigmaplus))';%% previous using 1, which is wrong
w_dS02=((RdS02pred)'/(Sigmaplus))';
w_dt=((Rdtpred)'/(Sigmaplus))';
w_gamma11=((Rgamma11)'/(Sigmaplus))';
w_gamma22=((Rgamma22)'/(Sigmaplus))';
w_gamma12=((Rgamma12)'/(Sigmaplus))';





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




