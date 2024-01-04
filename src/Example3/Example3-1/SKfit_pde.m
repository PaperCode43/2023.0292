function model = SKfit_pde(X, Y, B, Vhat, gammaP,lambda,Xtest,r,sigma,dS,Vmatrix)

[k, d] = size(X);
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

lbtau2 = 0.01;
%ubtau2 = 100;

% lbtheta = 0.001*ones(d,1); 

lbtheta = [0.0001;0.0001];
%ubtheta = 10*ones(d,1);

lb = [lbtau2;lbtheta];
% no upper bounds
ubtau2=200;
ubtheta=[10;10];
ub =[ubtau2;ubtheta];
%ub=[100;1;1;1];
%ub=[];


myopt = optimset('MaxFunEvals',100000000,'MaxIter',50000000);
% % % myopt = optimset('Display','iter','MaxFunEvals',1000000,'MaxIter',500,'GradObj','on');    %%% add gradient by SHH
%[in1,in2]=neg_LL_p([tau2_0;theta_0],k,d,D,B,V,Y,gammaP,lambda);
parms = fmincon(@(x) neg_LL_p_parms(x,k,d,D,B,V,Y,gammaP,lambda),...
        [tau2_0;theta_0],[],[],[],[],lb,ub,[],myopt); 


% record MLEs for tau2 and theta 
tau2hat = parms(1);
thetahat = parms(2:length(parms));


% the neg_LL for original problem
model.negLL = neg_LL([tau2hat; thetahat],k,d,D,B,V,Y,gammaP);



% 
% %%% record the -log likelihood
% model.negLL = neg_LL([tau2hat; thetahat],k,d,D,B,V,Y,gammaP);


% MLE of beta is known in closed form given values for tau2, theta
% and is computed below
% calculate estimates of the correlation and covariance matrices
Rhat = corrfun(thetahat,D,gammaP);
Sigmahat  = tau2hat*Rhat + V;
% calculate betahat
temp1 = B'/Sigmahat;
temp2 = temp1*B;
temp3 = temp1*Y;
betahat = temp2\temp3;

% issue warnings related to constraints 
warningtolerance = 0.001;
%theta_int=thetahat;
if min(abs(lbtheta - thetahat)) < warningtolerance
    warning('thetahat was very close to artificial lower bound');
    theta_int=[0.01;0.5];
    %theta_int=thetahat;
else
    theta_int=thetahat;
end
% if abs(lbtau2 - tau2hat) < warningtolerance
%     warning('tau2hat was very close to artificial lower bound');
% end
% 
if min(abs(ubtheta - thetahat)) < warningtolerance
    warning('thetahat was very close to artificial upper bound');
    theta_int=[0.01;0.2];
    %theta_int=thetahat;
end
% if abs(ubtau2 - tau2hat) < 100*warningtolerance
%     warning('tau2hat was very close to artificial upper bound');
%     theta_int=[0.01;0.01;0.01];
%     tau2hat=20;
% end

%% second stage optimization
%par_opt=fmincon(@(x) obj2(x,X,Xtest,Vhat,D,gammaP,a,Y,betahat),...
      %  [tau2hat;thetahat],[],[],[],[],0.1*lb,ub,[],myopt);
%if sum(Vhat)==0
%     par_opt=fmincon(@(x) obj2_nonoise1(x,X,Xtest,Vhat,D,gammaP,Y,betahat,tau2hat,r,sigma),...
%         thetahat,[],[],[],[],lbtheta,ub,[],myopt);
%     tau2opt=tau2hat;
%     thetaopt = par_opt;
%B1=polybasis(X,0);
%Ynew=Y - B1*betahat;
Yp=[Y;dS];





    
%     theta_int=thetahat/3;
%     m=10;J=100;T=1;
%     nopde_ob=obj2_whole_3D(thetahat,X,Xtest,Vhat,D,gammaP,Yp,betahat,tau2hat,r,sigma,Rho,T,m,J,Vmatrix);
%     par_optm=zeros(3,10);
%     pde_obm=10^18*ones(10,1);
%     for uu=1:6
%         %uu
%     %par_opt=fmincon(@(x) obj2_whole_3D(x,X,Xtest,Vhat,D,gammaP,Y,betahat,tau2hat,r,sigma,Rho,T,m,J),...
%         %theta_int,[],[],[],[],lbtheta,ub,[],myopt);
%     par_opt=fmincon(@(x) obj2_whole_3D(x,X,Xtest,Vhat,D,gammaP,Yp,betahat,tau2hat,r,sigma,Rho,T,m,J,Vmatrix),...
%         theta_int,[],[],[],[],[0.001;0.001;0.001],[1;1;1],[],myopt);
%     %%%wether theta is good
%     pde_ob=obj2_whole_3D(par_opt,X,Xtest,Vhat,D,gammaP,Yp,betahat,tau2hat,r,sigma,Rho,T,m,J,Vmatrix);
%     par_optm(:,uu)=par_opt;
%     pde_obm(uu)=pde_ob;
%     if pde_ob<0.8*nopde_ob
%         thetaopt = par_opt;
%         break
%     else
%         theta_int=theta_int+max([thetahat/3,2*lbtheta],[],2);
%         %myopt = optimset('MaxFunEvals',10000,'MaxIter',5000);
%         min_pde_ob=min(pde_obm);
%         idxx=find(pde_obm==min_pde_ob);
%         %thetaopt = par_optm(:,idxx);
%         thetaopt=thetahat; %% change 2022-11-21
%     end
%     end

    %theta_int=thetahat;
    Vmatrix1=Vmatrix+0.00001*diag(ones(size(Vmatrix,1),1));
    m=10;J=100;T=1;
    nopde_ob=obj2_whole_2D(thetahat,X,Xtest,Vhat,D,gammaP,Yp,betahat,tau2hat,r,sigma,T,m,J,Vmatrix1);
 
    
    myopt2 = optimset('MaxFunEvals',100000,'MaxIter',50000);
    par_opt=fmincon(@(x) obj2_whole_2D(x,X,Xtest,Vhat,D,gammaP,Yp,betahat,tau2hat,r,sigma,T,m,J,Vmatrix1),...
       theta_int,[],[],[],[],0.3*theta_int,3*theta_int,[],myopt2);
    %par_opt=fmincon(@(x) obj2_whole_2D(x,X,Xtest,Vhat,D,gammaP,Yp,betahat,tau2hat,r,sigma,T,m,J,Vmatrix1),...
       %3*theta_int,[],[],[],[],lbtheta,ubtheta,[],myopt2);



    
    
    %%%wether theta is good
    pde_ob=obj2_whole_2D(par_opt,X,Xtest,Vhat,D,gammaP,Yp,betahat,tau2hat,r,sigma,T,m,J,Vmatrix1);
   
    if pde_ob<nopde_ob
        thetaopt = par_opt;
    else
        thetaopt = thetahat; %% change 2023-11-26
    end
   
    
    
    tau2opt=tau2hat;
%     par_opt2d=fmincon(@(x) obj2_whole(x,X,Xtest,Vhat,D,gammaP,Y,betahat,tau2hat,r,sigma,T,m,J),...
%          thetahat,[],[],[],[],lbtheta,ub,[],myopt);
    
    %tau2opt=tau2hat;
%     par_opt2d=fmincon(@(x) obj2_whole(x,X,Xtest,Vhat,D,gammaP,Y,betahat,tau2hat,r,sigma,T,m,J),...
%         thetahat,[],[],[],[],lbtheta,ub,[],myopt)
%objstage2_integral=obj2_whole(par,X,Xtest,Vhat,D0,gammaP,Y,beta,tau2,r,sigma,T,m,J)

%else
%    par_opt=fmincon(@(x) obj2_noise(x,X,Xtest,Vhat,D,gammaP,b,Y,betahat),...
%        [tau2hat;thetahat],[],[],[],[],0.1*lb,ub,[],myopt);
%    tau2opt = par_opt(1);
%    thetaopt = par_opt(2:length(par_opt));
%end
model.negLL = neg_LL([tau2opt; thetaopt],k,d,D,B,V,Y,gammaP);
Rhatopt = corrfun(thetaopt,D,gammaP);
Sigmahat_opt = tau2opt*Rhatopt + V;
% calculate betahat
temp1 = B'/Sigmahat_opt;
temp2 = temp1*B;
temp3 = temp1*Y;
betahatopt = temp2\temp3;



%betahatopt=betahat;

%l0=obj2_nonoise1(lbtheta,X,Xtest,Vhat,D,gammaP,Y,betahatopt,tau2hat,r,sigma)
%l0=obj2_whole(lbtheta,X,Xtest,Vhat,D,gammaP,Y,betahatopt,tau2hat,r,sigma,T,m,J)

%l1=obj2_nonoise1(thetaopt,X,Xtest,Vhat,D,gammaP,Y,betahatopt,tau2hat,r,sigma)
%l1=obj2_whole(thetaopt,X,Xtest,Vhat,D,gammaP,Y,betahatopt,tau2hat,r,sigma,T,m,J)
%l2=obj2_nonoise1(thetahat,X,Xtest,Vhat,D,gammaP,Y,betahatopt,tau2hat,r,sigma)
%l2=obj2_whole(thetahat,X,Xtest,Vhat,D,gammaP,Y,betahatopt,tau2hat,r,sigma,T,m,J)
% output MLEs and other things useful in prediction
model.tausquared =  tau2opt;
model.beta = betahatopt;
%model.beta = betahat;%%2022-11-22
model.theta = thetaopt;
model.X = X;
% model.minX = minX;
% model.maxX = maxX;
model.gamma = gammaP;
model.Sigma = Sigmahat_opt;
model.Z = Y-B*betahat;
model.B = B;
model.parms=parms;


%%%%%% print parameters
fprintf('\n');
fprintf('betahat_pde = ');
if length(betahat)>1
    fprintf('[');
end
for i=1:length(betahat)
    fprintf('%.4f; ',betahatopt(i));
end
if length(betahat)>1
    fprintf(']');
end
fprintf('\n');
fprintf('tau2hat_pde = %.4f  \n',tau2opt);
fprintf('thetahat_pde = ');
if length(thetahat)>1
    fprintf('[');
end
for i=1:length(thetahat)
    fprintf('%.4f; ',thetaopt(i));
end
if length(thetahat)>1
    fprintf(']');
end
fprintf('\n');
fprintf('log-likelihood = %.4f  \n',-model.negLL);
fprintf('AIC = %.4f  \n',-2*(-model.negLL)+2*(1+length(betahat)+length(thetahat)));
fprintf('\n');
%fprintf('nopde_ob = %.4f  \n',nopde_ob);
%fprintf('uu = %.4f  \n',uu);
%fprintf('pde_ob = %.4f  \n',pde_ob);

