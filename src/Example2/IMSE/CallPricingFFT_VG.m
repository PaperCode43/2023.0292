 

function call_price_fft = CallPricingFFT_VG(n,S,K,T,r,d,sigma,nu,theta)
tic
lnS = log(S);
lnK = log(K);

%optAlpha = optimalAlpha(model,lnS,lnK,T,r,d,varargin{:});
optAlpha = .75;

DiscountFactor = exp(-r*T);


FFT_N = 2^n;                               % must be a power of two (2^14)
FFT_eta = 0.05;                             % spacing of psi integrand


FFT_lambda = (2 * pi) / (FFT_N * FFT_eta);  %spacing for log strike output (23)
FFT_b = (FFT_N * FFT_lambda) / 2;           % (20)

uvec = 1:FFT_N;
%log strike levels ranging from lnS-b to lnS+b
ku = - FFT_b + FFT_lambda * (uvec - 1);     %(19)

jvec = 1:FFT_N;
vj = (jvec-1) * FFT_eta;

tmp = DiscountFactor * psi(vj,optAlpha,lnS,T,r,d,sigma,nu,theta) .* exp(1i * vj * (FFT_b)) * FFT_eta;
tmp = (tmp / 3) .* (3 + (-1).^jvec - ((jvec - 1) == 0) );   %applying simpson's rule
cpvec = real(exp(-optAlpha .* ku) .* fft(tmp) / pi);        %call price vector resulting in equation 24

indexOfStrike = floor((lnK + FFT_b)/FFT_lambda + 1); 
iset = max(indexOfStrike)+1:-1:min(indexOfStrike)-1;
xp = ku(iset);
yp = cpvec(iset);
call_price_fft = real(interp1(xp,yp,lnK));
toc
end

%analytical formula for zhi in equation ( 6 ) of Madan's paper
function ret = psi(v,alpha,lnS,T,r,d,sigma,nu,theta)
% Variance Gamma
    u=v - (alpha + 1) * 1i;
    omega = (1/nu)*( log(1-theta*nu-sigma*sigma*nu/2) );
    tmp = 1 - 1i * theta * nu * u + 0.5 * sigma * sigma * u .* u * nu;
    y = 1i * u * (lnS + (r + omega - d) * T ) - T*log(tmp)/nu;
    ret=exp(y)./ (alpha.^2 + alpha - v.^2 + 1i * (2 * alpha + 1) .* v);
%end
  %ret = exp(feval(@CharacteristicFunctionLib, model, v - (alpha + 1) * 1i,varargin{:})) ./ (alpha.^2 + alpha - v.^2 + 1i * (2 * alpha + 1) .* v);
end

% function ret = psialpha(model,alpha,varargin)
%   ret = exp(feval(@CharacteristicFunctionLib, model, - (alpha + 1) * 1i,varargin{:}))./ (alpha.^2 + alpha);
% end