function F=SPncf(xords,n1,n2,theta1,theta2,acclevel)
% F=SPncf(xords,n1,n2,theta1,theta2=0,acclevel=2)
% Computes the saddlepoint cdf of the doubly noncentral F distribution
% of (vector) xords with n1 and n2 degrees of freedom and noncentrality params theta1,theta2 
% if acclevel=2 [default], then add the Daniels (1987) term
%
%        U1 / n1
%   F = ---------  , where Ui = Yi'Yi, Yi ~ N(mui, I(ni))
%        U2 / n2     and thetai = mui'mui , i=1,2
%
% Marc Paolella, 1999.


if nargin<6, acclevel=2; end
if nargin<5, theta2=0; end

epsilon=1e-4; % how close can s get to zero before we switch to the limiting value at the mean.
tau=1e7;      % how large can 1/u^3 and 1/w^3 be before we switch back to regular LugRic.

flip=0;
%if ((theta2>0) & (theta1==0))
%  flip=1;
%  xords=1./xords; theta1=theta2; theta2=0;
%  temp=n1; n1=n2; n2=temp;
%else
%  flip=0;  
%end

xlen=length(xords);
F=zeros(xlen,1);
for xloop=1:xlen, 
  x=xords(xloop);
  l1= n2/n1; l2=-x;
  if (theta2==0)
    x2=x^2; a=n1; a2=n1^2; a3=a*a2; b=n2; b2=n2^2; t=theta1; t2=t^2; % simpler
    y=x2 * a3 + 2*x2*a2*t + 2*a2*x*b + 4*x2*a*b*t + a*t2*x2 + 2*a*t*x*b + b2*a + 4*x*b2*t;
    num=x*a*(a+2*b+t) - a*b - sqrt(a*y); den=4*b*x*(a+b);
    s=num/den; roots=s;
  else
    x2=x^2; n12=n1^2; n22=n2^2; %n13=n1^3; 
    n23=n2^3; t1=theta1; t2=theta2;
    a=8*(n23*x2+n22*n1*x2);
    b=4*(-2*n22*n1*x2 + x*n23 + 2*x*n1*n22 + x*t2*n22 - n12*n2*x2 - n1*n2*t1*x2);
    c=2*(n22*n1 + n12*n2*x2) - 4*(x*n12*n2 + n1*n2*t1*x + x*t2*n1*n2 + x*n1*n22);
    d=-n12*n2 + x*n12*n2 + x*t2*n12 - n1*n2*t1;
    
    % use Abrahamowitz and Stegun notation
    a2=b/a; a1=c/a; a0=d/a;
    q=a1/3 -a2^2/9; r=(a1*a2-3*a0)/6 - a2^3/27;
    m=q^3 + r^2;
    if m >= 0, error ('this should not happen!'), end
    s1=(r+sqrt(m))^(1/3); s2=(r-sqrt(m))^(1/3);
    sps=s1+s2; sms=s1-s2;
    % z1=sps - a2/3; z2=-sps/2 - a2/3 + sqrt(-3)*sms/2;
    z3=-sps/2 - a2/3 - sqrt(-3)*sms/2;
    roots=z3;
    
  end
  s=roots;
  v1 = 1/(1-2*s*l1); v2 = 1/(1-2*s*l2);
  K=0.5*(n1*log(v1) + n2*log(v2)) + s*(l1*theta1*v1 + l2*theta2*v2);
  kpp=2*(l1^2*v1^2*(n1+2*theta1*v1)  +  l2^2*v2^2*(n2+2*theta2*v2) );
  
  badpt=(K>0) | (kpp<=0); 
  if (abs(s)< epsilon) || badpt
    
    % numerically integrate the spa pdf
    if 1==1
      disp(['NOTE: For x=',num2str(x),', switching method to numeric integration of the SPA to the pdf'])
      cdf = spncfnumint(x,n1,n2,theta1,theta2);
    end  

    % other, not-so-great ideas.
    if 1==2  % Take s=0  and x = E[X_f]. 
             % Thus, same cdf value returned for all x in neighborhood s.t. abs(s)< S_tol
      l2=-(1+theta1/n1) / (1+theta2/n2);
      Kpp0=2*(l1^2*(n1+2*theta1)  +  l2^2*(n2+2*theta2) );
      Kppp0=8*(l1^3*(n1+3*theta1)  +  l2^3*(n2+3*theta2) );
      cdf=0.5+Kppp0/6/sqrt(2*pi)/Kpp0^(3/2);
    elseif 1==2 % Take s=0 but x is as given by the user, via l2.
      Kpp0=2*(l1^2*(n1+2*theta1)  +  l2^2*(n2+2*theta2) );
      Kppp0=8*(l1^3*(n1+3*theta1)  +  l2^3*(n2+3*theta2) );
      cdf=0.5+Kppp0/6/sqrt(2*pi)/Kpp0^(3/2);
    elseif 1==2   % Now use computed shat instead of zero in the limiting LugRic equation.
      K3=  8*(l1^3*v1^3*(n1+3*theta1*v1)  +  l2^3*v2^3*(n2+3*theta2*v2) );
      cdf=0.5+K3/6/sqrt(2*pi)/kpp^(3/2);
    end  
      
  else
    w=sign(s)*sqrt(-2*K); u=s*sqrt(kpp);
    npdf=normpdf(w);
    if acclevel==2
      npdf2 = npdf;
      K3=  8*(l1^3*v1^3*(n1+3*theta1*v1)  +  l2^3*v2^3*(n2+3*theta2*v2) );
      K4= 48*(l1^4*v1^4*(n1+4*theta1*v1)  +  l2^4*v2^4*(n2+4*theta2*v2) );
      K5=384*(l1^5*v1^5*(n1+5*theta1*v1)  +  l2^5*v2^5*(n2+5*theta2*v2) );
      kap3= K3/(kpp)^(3/2); 
      kap4= K4/(kpp)^(4/2); 
      kap5= K5/(kpp)^(5/2);
      
      tt1 = (kap4/8-5*kap3^2/24)/u;
      tt2 = 1/u^3;       bad=abs(tt2) > tau;
      tt3 = kap3/2/u^2;  
      tt4 = 1/w^3;       bad=(bad | (abs(tt4) > tau));
      term = tt1 - tt2 - tt3 + tt4;
      if bad
        if 1==2 % do nothing, keep the 2nd order term
          %term=term;
        elseif 1==1 % just drop the 2nd order term and stick with usual Lug-Ric formula
          term=0; 
        else    % try the other formula in Daniels
          %npdf2 = 1/sqrt(2*pi); % the limiting value
          term = kap5/40 - (5/48) * kap3*kap4 + (35/432) * kap3^2;
        end  
      end;
      
      O1=-npdf2 * term;
    else O1=0; 
    end
    cdf=normcdf(w)+npdf*(1/w - 1/u) + O1;
  end
  F(xloop)=cdf;
end  

if flip==1, F=1-F; end
F=F';

function F=spncfnumint(x,n1,n2,theta1,theta2)
% F=spncfnumint(x,n1,n2,theta1,theta2)
% cdf of the doubly noncentral F via numeric integration of the SPA pdf
% takes only scalar values of x

area=quadl(@SPncfpdf,0,x, 1e-8, 0 ,n1, n2, theta1, theta2);

up=2*x;
change=quadl(@SPncfpdf,x,up, 1e-8, 0 ,n1, n2, theta1, theta2);
new=area+change;
while change>1e-6
  lo=up;
  up=1.5*up;
  change=quadl(@SPncfpdf,lo,up, 1e-8, 0 ,n1, n2, theta1, theta2);
  new=new+change;
end
F=area/new;

