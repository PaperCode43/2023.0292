function model = SKfit_der(X, Y, B, Vhat, gammaP,lambda)
% fit a stochastic kriging model to simulation output
% k - number of design points
% d - dimension of design variables
% X - design points (k x d) 
% Y - (k x 1) vector of simulation output at each design point
% B - (k x b) matrix of basis functions at each design point
%     The first column must be a column of ones!
% Vhat - (k x 1) vector of simulation output variances
% Types of correlation function used to fit surface:
%    gammaP = 1: exponential
%    gammaP = 2: gauss
%    gammaP = 3: cubic
% 
% Examples
%       skriging_model = SKfit(X,Y,ones(k,1),Vhat,2);
% Use SK model with gauss correlation function to fit data, (X,Y,Vhat)
% X is design points, Y and Vhat are outputs at design points 
% (Y is mean, Vhat is intrinsic variance), non-trend estimate B = ones(k,1)

[k, d] = size(X);

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
D = zeros(k,k,d);
for p=1:d
    temp = zeros(k);
    if (gammaP == 3)
        temp(IdistX) = abs(distX(:,p));
    else
        temp(IdistX) = -abs(distX(:,p)).^gammaP;
    end
    D(:,:,p) = temp+temp';
end

% diagonal intrinsic variance matrix
V = diag(Vhat);
scaling = 1;



% initialize parameters theta, tau2 and beta
% inital extrinsic variance = variance of ordinary regression residuals
% if B is multi-variable, use lasso to overcome uninversible %% changed by Jiang

[~,bd]=size(B);
if bd ==1
betahat = (B'*B)\(B'*Y);
else
B_new=B(:,2:end);
[betahatB,fitinfo]=lasso(B_new,Y,'CV',5,'Lambda',0);
betahatC=fitinfo.Intercept;
betahat=[betahatC;betahatB];
end
tau2_0 = var(Y-B*betahat);
% see stochastic kriging tutorial for explanation
% make correlation = 1/2 at average distance
average_distance = mean(abs(distX));
theta_0 = zeros(d,1);
if gammaP == 3
    theta_0(1:d) = 0.5;
else
    theta_0(1:d) = (log(2)/d)*(average_distance.^(-gammaP));
end

% lower bounds for parameters tau2, theta
% naturally 0; increase to avoid numerical trouble

% lbtau2 = 0.01;
% 
% % lbtheta = 0.001*ones(d,1); 
% 
% lbtheta = 0.01*ones(d,1);
% 
% lb = [lbtau2;lbtheta];
% % no upper bounds
% ub =[];

lbtau2 = 0.1;
ubtau2 = 100;

%lbtheta = 0.01*ones(d,1); 

lbtheta = [0.02;0.2;0.2;0.1];
%ubtheta = 10*ones(d,1);
ubtheta = [2;50;50;10];

lb = [lbtau2;lbtheta];
% no upper bounds
ub =[ubtau2;ubtheta];


% maximize profile log-likelihood function (-"neg_LL"), i.e.,
% minimize the - log-likelihood funtion ("neg_LL"),
% subject to lower bounds on tau2 and theta
% GA is more suitable for high dimension problem


%%%% GA algorithm (genetic algorithm), requires function with row vector as input !
% for reproducibility
seed = 100; 
s_pseudo = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(s_pseudo);

% options = gaoptimset('Display','iter');

 %options = gaoptimset('Display','iter', 'Generations',200);
options = gaoptimset('Generations',200);

%warning off
parms = ga(@(x) neg_LL(x,k,d,D,B,V,Y,gammaP), 1+d, [],[],[],[],lb,ub,[],[],options);
parms = parms';  


% use the solution  from GA as initial solution
tau2_0 = parms(1);
theta_0 = parms(2:length(parms));

%warning on

%%%% fmincon algorithm
% myopt = optimset('Display','iter','MaxFunEvals',1000000,'MaxIter',500);
%myopt = optimset('MaxFunEvals',1000000,'MaxIter',5000,'GradObj','on');
myopt = optimset('MaxFunEvals',10000000,'MaxIter',50000);
% % % myopt = optimset('Display','iter','MaxFunEvals',1000000,'MaxIter',500,'GradObj','on');    %%% add gradient by SHH
%[in1,in2]=neg_LL_p([tau2_0;theta_0],k,d,D,B,V,Y,gammaP,lambda);
parms = fmincon(@(x) neg_LL_p_parms(x,k,d,D,B,V,Y,gammaP,lambda),...
        [tau2_0;theta_0],[],[],[],[],lb,ub,[],myopt); 


% record MLEs for tau2 and theta 
tau2hat = parms(1);
thetahat = parms(2:length(parms));


% the neg_LL for original problem
model.negLL = neg_LL([tau2hat; thetahat],k,d,D,B,V,Y,gammaP) - k*log(scaling);


%%%% scale back to original parameter
if scaling ~= 1
    tau2hat = tau2hat/(scaling^2);
    Y = Y/scaling;
    B = B/scaling;
    V = V/(scaling^2);
end


% 
% %%% record the -log likelihood
% model.negLL = neg_LL([tau2hat; thetahat],k,d,D,B,V,Y,gammaP);


% MLE of beta is known in closed form given values for tau2, theta
% and is computed below
% calculate estimates of the correlation and covariance matrices
Rhat = corrfun_der(thetahat,-D,gammaP,1);%% partial_idx=1 becasue X is 1-D
Sigmahat  = tau2hat*Rhat + V;
% calculate betahat
temp1 = B'/Sigmahat;
temp2 = temp1*B;
temp3 = temp1*Y;
betahat = temp2\temp3;



% issue warnings related to constraints 
warningtolerance = 0.001;
if min(abs(lbtheta - thetahat)) < warningtolerance
    warning('thetahat was very close to artificial lower bound');
end
if abs(lbtau2 - tau2hat) < warningtolerance
    warning('tau2hat was very close to artificial lower bound');
end

% output MLEs and other things useful in prediction
model.tausquared =  tau2hat;
model.beta = betahat;
model.theta = thetahat;
model.X = X;
% model.minX = minX;
% model.maxX = maxX;
model.gamma = gammaP;
model.Sigma = Sigmahat;
model.Z = Y-B*betahat;
model.B = B;


%%%%%% print parameters
fprintf('\n');
fprintf('betahat = ');
if length(betahat)>1
    fprintf('[');
end
for i=1:length(betahat)
    fprintf('%.4f; ',betahat(i));
end
if length(betahat)>1
    fprintf(']');
end
fprintf('\n');
fprintf('tau2hat = %.4f  \n',tau2hat);
fprintf('thetahat = ');
if length(thetahat)>1
    fprintf('[');
end
for i=1:length(thetahat)
    fprintf('%.4f; ',thetahat(i));
end
if length(thetahat)>1
    fprintf(']');
end
fprintf('\n');
fprintf('log-likelihood = %.4f  \n',-model.negLL);
fprintf('AIC = %.4f  \n',-2*(-model.negLL)+2*(1+length(betahat)+length(thetahat)));


see = [-2*(-model.negLL)+2*(1+length(betahat)+length(thetahat)), tau2hat];


%%%%%% hypothesis test statistic
if length(betahat) == 2      %%% hypothesis test _ add by SHH
    m22 = (B'/Sigmahat*B)^(-1);
    Z = betahat(2)/sqrt(m22(2,2));    
    pV = 2* normcdf(-abs(Z));
    fprintf('Usefullness Test: Z = %.4f,  P_Value = %.4f  \n',Z,pV);
    
    see(1,3) = pV;
end

model.see = see;


