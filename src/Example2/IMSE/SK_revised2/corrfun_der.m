function R = corrfun_der(theta,D,gammaP,partial_idx)

% theta = scale parameter, (d*1)
% D = distance matrix (k1*k2*d)
% R = cubic correlation function (k1*k2)

if gammaP ~= 2   
    % return cubic correlation function
    error('This is correlation of derivative surface, gamma=2!')   
else
    % return gauss or exponential correlation function
    d = size(theta,1);
    k1 = size(D,1);
    k2 = size(D,2);
    Rl = exp(sum(-D.*repmat(reshape(theta,[1 1 d]),[k1 k2]),3));
    pR = D(:,:,partial_idx);
    R=2 * theta(partial_idx) * (1 - 2 * theta(partial_idx) * pR).*Rl;
end