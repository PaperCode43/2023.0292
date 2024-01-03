function [f,neg_gradient] = neg_LL(parms,k,d,D,B,V,Y,gammaP)

%%%% convert to column vector (since GA use row vector as input!)
parms = parms(:);

% negative of log profile likelihood  
% This is a function of theta and tau2.
% we have profiled out over beta
% parms = [tausq;theta] 
% k = number of design points
% d = dimension of space
% Dist = matrix of distances between design points
% B = basis functions
% V = intrinsic covariance matrix
% Y = simulation output

%%%% conditon check
if ((length(parms)~=d+1))
    error('Parameter vector must be length d+1');
end
if not(all(size(Y)==[k 1]))
    error('Output vector must be k by 1.')
end
if (size(B,1)~=k)
    error('Basis function matrix must have k rows.')
end
if (not(all(size(V)==[k k])))
    error('The intrinsic covariance matrix must be k by k.')
end

if(parms(1) < 0.001 || min(parms(2:d+1)) < 0.001)
    f = inf;
    fprintf('The parameters reach the boundary \n');
    return;
end


tausq = parms(1);
theta = parms(2:d+1);

% get correlation matrix given theta
R = corrfun(theta,D,gammaP);
 
% sum of extrinsic and intrinsic covariances
Sigma  = tausq*R + V;

% % the optimal beta given theta and tau2
% temp1 = B'/Sigma;
% temp2 = temp1*B;
% temp3 = temp1*Y;
% beta = temp2\temp3;
% 
% % log likelihood
% f = 0.5*k*log(2*pi) + 0.5*log(det(Sigma)) + 0.5*(Y-B*beta)'/Sigma*(Y-B*beta);


% more numerically stable
[U,pd] = chol(Sigma);
if(pd>0)
    error('covariance matrix is nearly singular');
end

% invert it via Cholesky factorization
L = U';
Linv = inv(L);
Sinv = Linv'*Linv;

% the optimal beta given theta and tau2

beta = (B'*Sinv*B)\(B'*(Sinv*Y));
Z = L\(Y-B*beta);


% negative log likelihood
f = 0.5*k*log(2*pi) + log(det(L)) + 0.5*Z'*Z;
% Note:  log(det(L)) ==  0.5* log(det(Sigma))


% gradient required, only defined when theta is 1-D
if nargout > 1 
    
    A = Sigma;
    F = B;
    C = R;
    C1 = sum(D,3).*C;   %%%%  D has already incorporated a '-'
    tau2 = tausq;     
        
    beta_par_tau2 = (F'*inv(A)*F) \ F'*(inv(A)*C*inv(A)) * (F/(F'*inv(A)*F)*F'*inv(A)*Y - Y);    
    beta_par_theta = tau2* ((F'*inv(A)*F) \ F'*(inv(A)*C1*inv(A))) * (F/(F'*inv(A)*F)*F'*inv(A)*Y - Y);

    % this is the gradient of log-likelihood function w.r.t. tau^2, theta_1, ..., theta_d
    gradient_LL = ...
      [-0.5*trace(inv(A)*C) + 0.5*(Y-F*beta)' *inv(A) * C * inv(A)*(Y-F*beta) + (Y-F*beta)'*inv(A)*F * beta_par_tau2; ...
       -0.5*trace(tau2*inv(A)*C1) + 0.5*tau2*(Y-F*beta)' *inv(A)*C1*inv(A)*(Y-F*beta) + (Y-F*beta)'*inv(A)*F * beta_par_theta];
    neg_gradient = - gradient_LL;
end

end
