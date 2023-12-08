function [mu, Sigma] = DCC1step(data)
% INPUT: data is a T X d matrix of asset log returns
% Output: 
%   mu: the mean of returns, 1 x d vector
%   sigma: one step ahead prediction of the asset covariance matrix
%     computed using the Gaussian DCC-GARCH model

profile=10; [T,d]=size(data); 
mu=mean(data);  data_demean=data-repmat(mu,T,1);
garchP=zeros(3,d); sigmamat = zeros(T,d); GARCHresid = zeros(T,d);
for n = 1:d
  garchP(:,n) = normalGARCHPL(data_demean(:,n),[],profile);
  [GARCHresid(:,n),sigmamat(:,n)]=ungarch(data_demean(:,n),garchP(:,n));
end
%DCC: Q_t=S*(1-a-b)+a*eps_{t-1}eps'_{t-1} +b*Q_{t-1}
  initvec=[]; S=[]; Gamma0=[]; 
  [a,b,Smat,~,Gamma0] = DCCestimate(GARCHresid,initvec,S,Gamma0);
  Rmat  = zeros(d,d,T); S=Smat; Q=Gamma0;
  for t=1:T
    if t>1
      Q = DCCengine(GARCHresid(t-1,:),a,b,S,Q);
      Corrmat = diag(sqrt(diag(Q)).^(-1))*Q*diag(sqrt(diag(Q)).^(-1));
      Corrmat = SymPDcorrmatrix(Corrmat);
    else
      Corrmat = S; Q=S; Corrmat=SymPDcorrmatrix(Corrmat);
    end
    Rmat(:,:,t)=Corrmat;
  end
  Q_f=(1-a-b)*Rmat(:,:,1)+a*GARCHresid(end,:)'*GARCHresid(end,:)+b*Rmat(:,:,end);
  R_f=diag(diag(Q_f).^(-0.5))*Q_f*diag(diag(Q_f).^(-0.5));
  H_f=diag(sigmamat(end,:))*Rmat(:,:,end)*diag(sigmamat(end,:));
% Recombine GARCH11 and Correlation matrix forecast
D_f=diag(sqrt(garchP(1,:)+garchP(2,:).*GARCHresid(end,:).^2+garchP(3,:).*diag(H_f)')); 
Sigma=D_f*R_f*D_f;  Sigma=SymPDcovmatrix(Sigma);


function [a,b,Smat,ll,Gamma0] = DCCestimate(resid,initvec,S,Gamma0)
fmintol=1e-6; MaxFunEvals=5000;
    opt = optimset('Display','off','TolFun',fmintol, ...
    'TolX',fmintol,'MaxFunEvals',MaxFunEvals,'LargeScale','Off');
%             C            D
bound.lo=    [0.001        0         ];
bound.hi=    [1-(1e-04)    1-(1e-04) ];
bound.which= [1            1         ];
if isempty(initvec), initvec=[0.05 0.93]; end
if isempty(S), S=corr(resid); end;
S=SymPDcorrmatrix(S); if isempty(Gamma0), Gamma0=S; end;
Gamma0 = SymPDcorrmatrix(Gamma0);
[pout,~] = fminunc(@(param) loglikDCC(param,resid,S,Gamma0,bound), ...
    einschrk(ab2CD(initvec),bound),opt);
param=CD2ab(einschrk(real(pout),bound,1)); a=param(1); b=param(2); Smat=S;
ll= -1 * loglikDCC(ab2CD(param),resid,S,Gamma0,[]);

function loglik=loglikDCC(param,resid,S,Q,bound)
if ~isempty(bound), param=einschrk(real(param),bound,999); end
a = param(1) * param(2); b = param(2) - param(1) * param(2);
[T,d] = size(resid); ll=zeros(T,1);
for t=1:T
  if t>1
   Q = DCCengine(resid(t-1,:),a,b,S,Q);
   Corrmat = diag(sqrt(diag(Q)).^(-1))*Q*diag(sqrt(diag(Q)).^(-1));
   Corrmat = SymPDcorrmatrix(Corrmat);
  else
   Corrmat=S; Q=S; Corrmat=SymPDcorrmatrix(Corrmat);
  end
  detGamma=det(Corrmat); Gammainv=Corrmat\eye(d);
  ll(t) = - 0.5 *(log(abs(detGamma)) + resid(t,:) * Gammainv * resid(t,:)');
end
loglik= - sum(ll);

function Q=DCCengine(resid_t,a,b,S,Qin)
Mmat = resid_t'*resid_t; Mmat=SymPDcovmatrix(Mmat);
Q=(1-a-b) * S + a * Mmat + b * Qin;

function [eout,sigvec]=ungarch(epsi,garchin)
global delta
int=garchin(1); qterms=garchin(2); pterms=garchin(3);
e=(abs(epsi)).^delta; 
lambda=( ( 2^(delta/2) ) / sqrt(pi) ) * gamma((delta+1)/2);

% sinit =  E[sigvec_0^delta], einit = E|e_0^delta| = lambda * sinit.
sinit=mean((abs(epsi)).^delta); einit=lambda*sinit;

% do the recursion in sigvec^delta
sigvec=zeros(length(e),1); sigvec(1)=int+qterms*einit+pterms*sinit;
for i=2:length(e), sigvec(i)=int + qterms*e(i-1) + pterms*sigvec(i-1); end
if any(sigvec<=0), error('hello'), end
sigvec=sigvec.^(1/delta); eout=epsi./sigvec;

function A=SymPDcovmatrix(A,tol)
[n,m]=size(A);
if ~(n==m), error('Input matrix has to be a square matrix '), end
if nargin<2, tol=1e-04; end
A=(A+A')/2;
try
  [V,D]=eig(A); seig=diag(D); bad=find(seig<tol);
catch %#ok<CTCH>
  %bp=1; 
end
if ~isempty(bad), seig(bad)=tol; D=diag(seig); A=V*D*V';end

function A=SymPDcorrmatrix(A,tol)
[n,m]=size(A);
if ~(n==m), error('Input matrix has to be a square matrix '), end
if nargin<2, tol=1e-04; end
if sum(any(abs((A)-eye(n))>=1+(1e-02)))>0, warning('This is not a correlation matrix'); end
numCol = find( any( abs((A)-eye(n)) >= ( 1-(1e-16) ) ) );
numRow = find( any( abs((A)-eye(n))' >= ( 1-(1e-16) ) ) );
if ~isempty(numCol) || ~isempty(numRow)
  A(numRow,numCol) = sign(A(numRow,numCol)) * ( 1-(1e-16) ); 
  % the off-diagonal entries in (-1 1)
  warning('Some of the correlations were corrected');
end
A = A - diag(diag(A)) + eye(n,m); % ones on the diagonal
A = ( A + A' )/2;                 % symmetric matrix
[V,D] = eig(A); seig = diag(D); bad = find(seig < tol); % PD
if ~isempty(bad), seig(bad) = tol; D = diag(seig); A = V * D * V';end

function [ab,Vout] = CD2ab(CD,Vin)
%a=C*D;b=D-CD so 0<a<1, 0<b<1 and 0<a+b<1 can be achived by 0<C<1 and 0<D<1
[n1, n2]=size(CD);
if n1>2 || n2>2, warning('CD2ab programmed wrong'), end
ab(1) = CD(1).*CD(2); ab(2) = CD(2)-CD(1).*CD(2);
if nargout>1  % now adjust the standard errors
  J=[CD(2), CD(1); -CD(2), (1-CD(1))]; Vout=J * Vin *J;
end

function [CD, Vout] =ab2CD(ab,Vin)
%C=a/(a+b);D=a+b so 0<a<1, 0<b<1 and 0<a+b<1 can be achived by 0<C<1 and 0<D<1
[n1, n2]=size(ab);
if n1>2 || n2>2, warning('ab2CD programed wrong'), end %#ok<*WNTAG>
CD(1) = ab(1) ./ (ab(1) + ab(2)); CD(2) = ab(1) + ab(2);
if nargout>1  % now adjust the standard errors
  J=[1./(ab(1) + ab(2)) - ab(1)./(ab(1) + ab(2)).^2, -ab(1)/(ab(1) + ab(2)).^2 ; 1,1];
  Vout=J * Vin * J;
end

