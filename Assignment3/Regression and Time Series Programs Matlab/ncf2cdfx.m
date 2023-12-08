function ncf2cdfx=ncf2cdfx(alpha,n1,n2,theta1,theta2)
% finds the cutoff value of the (possibly doubly noncentral) F distribution
%   using the SPA.
% Compare to Matlab's built in ncfinv and finv.

show=0; % toggle between 1 and 0 to show (and not show) the iterative output of fmin

if (theta1>0) && (theta2>0) 
  xval = 1.5*theta1/theta2;
else
  xval=1;
end    

multip=1; cdf=2;
while (cdf>alpha)
  versuch= xval/multip;
  cdf = spncf(versuch,n1,n2,theta1,theta2);
  multip= multip*2;
end
lob= versuch;

multip= 1; cdf=-1;
while (cdf<alpha)
  versuch= xval*multip;
  cdf= spncf(versuch,n1,n2,theta1,theta2);
  multip= multip*2;
end
hib= versuch;

if 1==1
  ncf2cdfx=fminbnd(@(x) spncf_(x,n1,n2,theta1,theta2,alpha),lob,hib,optimset('TolX',1e-5,'Display','off'));
else % use bisection:
  versuch = (lob+hib)/2; valid=0; TOL=1e-8;
  while (valid~=1)
    cdf= spncf(versuch,n1,n2,theta1,theta2);
    valid= (abs(cdf-alpha)<TOL);
    if (valid==1), ncf2cdfx= versuch;
    else
      if (cdf<alpha), lob= versuch; else hib= versuch; end
      versuch= (lob+hib)/2;
    end
  end
end

end

function disc=spncf_(x,n1,n2,theta1,theta2,alpha)
disc=abs(spncf(x,n1,n2,theta1,theta2,2) - alpha);
end