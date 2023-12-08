function [f,F]=weightedsumofchisquare(xvec,a,n,theta)
% [f,F]=weightedsumofchisquare(xvec,a,n,theta)
% use inversion formulae for pdf and cdf of
% X=sum a_i X_i where X_i ~ chi2(n_i,theta_i)
% See Paolella (2007, page 351, eq. 10.19).

if nargin<4, theta=zeros(length(a),1); end
bordertol=1e-8; lo=bordertol; hi=1-bordertol; tol=1e-7;
xl=length(xvec); F=zeros(xl,1); f=F;
for loop=1:length(xvec)
  x=xvec(loop);
  f(loop)= quadl(@fff,lo,hi,tol,[],x,a,n,theta,1) / pi;
  if nargout>1
    F(loop)=0.5-(1/pi)* quadl(@fff,lo,hi,tol,[],x,a,n,theta,0);
  end
end;

function I=fff(uvec,x,a,n,theta,dopdf);
for ii=1:length(uvec)
  u=uvec(ii);  t  = (1-u)/u;
  v = 1./ (1-2*a*i*t);
  cf = exp( sum( n.*log(v)/2  + i*t*a.*theta.*v) );  % HERE IS THE CHARACTERISTIC FUNCTION
  z  = exp(-i*t*x) .* cf;
  if dopdf==1, g=real(z); else g=imag(z)./t; end
  I(ii)  = g / u^2;
end


