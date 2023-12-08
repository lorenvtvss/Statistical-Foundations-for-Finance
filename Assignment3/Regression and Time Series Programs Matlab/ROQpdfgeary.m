function [f]=ROQpdfgeary(rvec,A,B,mu,tol)
%f=ROQpdfgeary(r,A,B,mu,tol=1e-7)
%cdf and pdf of x'Ax/x'Bx, B>=0, x~N(mu,I), at the values in rvec.
if nargin<5,tol=1e-7;end
T=length(A);F=rvec*0;f=F;
if nargin<4,mu=zeros(T,1);end
A=.5*(A+A');B=.5*(B+B');

for loop=1:length(rvec)
    r=rvec(loop);
    C=A-r*B;
    [P,L]=eig(C);
    nu=P'*mu;
    wgt=diag(L);
    H=P'*B*P;
    if 1==1 % use the new matlab quadgk routine; it can handle singularities at the endpoints
      f(loop)=(1/pi) * quadgk(@(itvec)integrand(itvec,wgt,nu,H),0,1);
    else
      f(loop)=(1/pi) * quadl(@integrand,0,1,tol,[],wgt,nu,H);
    end
end

function [M] = integrand(itvec,wgt,nu,H)

M=zeros(size(itvec));
for loop=1:length(itvec)
    if itvec(loop)==0,
        itvec(loop)=0;
    else
        u=1./itvec(loop)-1;
        if 1==2,%do complex arithmetic
            t=i*u;
            iD=inv(eye(length(wgt))-2*t*diag(wgt));
            v=1./(1-2*t*wgt);
            G=trace(iD*H)+nu'*iD*H*iD*nu;
            M(loop)=real(exp(.5*sum(log(v))+t*sum(wgt.*nu.^2.*v)).*G)./itvec(loop).^2;
        else%imhof style real arithmetic
            nc=nu.^2;
            p = u*wgt; b=p.^2; d=1./(1+b);c=1+b;
            beta=sum(atan(p) + nc.*p.*d)/2;
            ss1=sum(nc.*b.*d);  ss2=-sum(log(d)); gam = exp(0.5*ss1 + 0.25*ss2);
            I=eye(length(wgt));
            L=diag(wgt);
            iF=diag(d);
            rho=trace(H*iF)+nu'*iF*(H-u^2*L*H*L)*iF*nu;
            delta=trace(H*L*iF)+2*nu'*iF*H*L*iF*nu;
            M(loop)=(cos(beta)*rho-u*sin(beta)*delta)/2/gam/itvec(loop)^2;
        end
    end
end