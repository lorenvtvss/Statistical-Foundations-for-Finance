function w = PortMNS(data, tau, DCC)
% Portfolio optimization using either:
%   the iid Markowitz setting 
%     or the DCC Markowitz setting
%   with non-negative portfolio weights, i.e., no short selling.
% INPUT:
% data is a T X d matrix of asset returns
% tau is the desired return, expressed as a yearly percentage 
%   (based on 250 trading days), e.g., tau=10.
% DCC is boolean, to use DCC or just the iid case (default).
% OUTPUT
%  w is the vector portfolio weight that satisfies the mean constraint 
%    and minimizes the variance.

if nargin<3, DCC=0; end
if DCC, [mu,Sigma]=DCC1step(data); else mu = mean(data); Sigma = cov(data); end
DEDR=100*((tau/100 + 1)^(1/250)-1); feas=max(mu) >= DEDR;
if feas, w=meanvar(mu,Sigma,DEDR)'; else w=zeros(length(mu),1); end

function w = meanvar(mu, Sigma, tau)
opt=optimset('Algorithm','active-set','LargeScale', 'off','Display','off'); 
d=length(mu); A = -mu; B = -tau; LB = zeros(1,d); UB = ones(1,d); w0=UB/d;
Aeq = ones(1,d); Beq = 1; % sum(w) = 1
w = fmincon(@fv, w0, A, B,  Aeq, Beq, LB, UB, [], opt, Sigma);

function f = fv(w, Sigma), f = w * Sigma * w';
