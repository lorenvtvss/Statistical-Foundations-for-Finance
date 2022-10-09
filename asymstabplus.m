% pdf and, optionally, cdf of the asymmetric stable. See also asymstabpdf.m

function [f,F] = asymstabplus(xvec,a,b,c,d)

if nargin<3, b=0; end
bordertol=1e-8; lo=bordertol ; hi=1-bordertol ; tol =1e-7;
xl=length(xvec); F=zeros(xl,1); f=F;
for loop=1:length(xvec)
    x=(xvec(loop)-d)/c;
    f(loop)=(integral(@(uvec) fffplus(uvec,x,a,b,c,d,1), lo, hi)/pi)/c;
    if nargout>1
        F(loop)=0.5-(1/pi)*integral(@(uvec) fffplus(uvec,x,a,b,c,d,0), lo, hi);
    end
end

