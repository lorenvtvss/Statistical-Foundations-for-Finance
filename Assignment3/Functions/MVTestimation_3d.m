function [param,stderr,iters,loglik,Varcov] = MVTestimation_3d(x, weights)
% param: (k, mu1, mu2, mu3, Sigma_11, Sigma_12, Sigma_13, Sigma_22, Sigma_23, Sigma_33)
[nobs d]=size(x); if d~=3, error('not done yet, use EM'), end
if d==3
    %%%%%%%%        k       mu1     mu2    mu3     s11     s12  s13   s22   s23     s33  
    bound.lo= [     0.2     -1      -1     -1     0.01    -90   0.01  -90  0.01    -90];
    bound.hi= [     20      1        1      1      90      90   90    90    90     90];
    bound.which=[   1       0        0      0      1       1    1     1     1       1];
    initvec =[      2      -0.8     -0.2   -0.8    20      2    10    20    2      10]; 
end

maxiter=300; tol=1e-7; MaxFunEvals=length(initvec)*maxiter;
opts=optimset('Display','iter','Maxiter',maxiter,'TolFun',tol,'TolX',tol,...
    'MaxFunEvals',MaxFunEvals,'LargeScale','Off');
[pout,fval,~,theoutput,~,hess]= ...
    fminunc(@(param) MVTloglik(param,x,bound, weights),einschrk(initvec,bound),opts);
V=inv(hess)/nobs; % Don't negate because we work with the negative of the loglik
[param,V]=einschrk(pout,bound,V); % transform and apply delta method to get V
param=param'; Varcov=V; stderr=sqrt(diag(V)); % Approximate standard errors
loglik=-fval*nobs; iters=theoutput.iterations;
end

function ll=MVTloglik(param,x,bound,weights)
if nargin<3, bound=0; end
if isstruct(bound), param=einschrk(real(param),bound,999); end
[nobs, d]=size(x); Sig=zeros(d,d); k=param(1); mu=param(2:4); % Assume d=2
Sig(1,1)=param(5); Sig(1,2)=param(6); Sig(1,3)=param(7); Sig(2,2)=param(8);Sig(2,3)=param(9); Sig(3,3)=param(10); Sig(2,1)=Sig(1,2); Sig(3,1)=Sig(1,3); Sig(3,2)=Sig(2,3); 
if min(eig(Sig))<1e-10, ll=1e5;
else
    pdf=zeros(nobs,1);
    for i=1:nobs, pdf(i) = mvtpdfmine(x(i,:),k,mu,Sig); end
    llvec=log(pdf); ll=-sum(llvec.*weights)/sum(weights); if isinf(ll), ll=1e5; end
    
end


function y = mvtpdfmine(x,df,mu,Sigma)
% x is a d X 1 vector. Unlike Matlab's version, cannot pass a matrix.
% Matlab's routine accepts correlation (not dispersion) matrix.
% So, just need to do the usual scale transfrom. For example:
%   x = [0.2 0.3]'; C = [1 .4; .4 1]; df = 2;
%   scalevec = [1 2]'; xx = x./scalevec; mvtpdf(xx,C,df)/prod(scalevec)
% Same as:
%   Sigma = diag(scalevec) * C * diag(scalevec); mvtpdfmine(x,df,[],Sigma)
d=length(x);
if nargin<3, mu = []; end, if isempty(mu), mu = zeros(d,1); end
if nargin<4, Sigma = eye(d); end
x = reshape(x,d,1); mu = reshape(mu,d,1); term = (x-mu)' * inv(Sigma) * (x-mu);
logN=-((df+d)/2)*log(1+term/df); logD=0.5*log(det(Sigma))+(d/2)*log(df*pi);
y = exp(gammaln((df+d)/2) - gammaln(df/2) + logN -logD);
end
end
