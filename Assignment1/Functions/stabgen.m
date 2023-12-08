%%% Generates a random sample of size nobs from S_(a,b)(d,c) distribution

function x=stabgen(nobs,a,b,c,d,seed)
if nargin<3, b=0; end, if nargin<4, c=1; end
if nargin<5, d=0; end, if nargin<6, seed=rand; end
z=nobs;
rand('twister',seed), V=unifrnd(-pi/2,pi/2,1,z);
rand('twister',seed+42), W=exprnd(1,1,z);
if a==1
    x=(2/pi)*(((pi/2)+b*V).*tan(V)-b*log((W.*cos(V))./((pi/2)+b*V)));
    x=c*x+d-(2/pi)*d*log(d)*c*b;
else
    Cab=atan(b*tan(pi*a/2))/(a);
    Sab=(1+b^2*(tan((pi*a)/2))^2)^(1/(2*a));
    A=(sin(a*(V+Cab)))./((cos(V)).^(1/a));
    B0=(cos(V-a*(V+Cab)))./W;
    B=(abs(B0)).^((1-a)/a);
    x=Sab*A.*(B.*sign(B0));
    x=c*x+d;
end

