function f = asymstabpdf(xvec,a,b,plotintegrand)
% pdf of the asymmetric stable . See also asymstab .m
% Set plotintegrand to 1, and xvec a scalar, to plot the integrand.

if nargin<4 , plotintegrand =0; end
if a==1, error('not ready yet') , end
xl = length(xvec); f = zeros(xl, 1); n = length(xvec);
for loop = 1:n
    x = xvec(loop);
    if (x < 0) , f(loop) = stab(-x , a, -b ) ; else f(loop) = stab(x,a,b) ; end
end
if plotintegrand && xl==1 % show a plot of the integrand
    if x<0 , x=-x ; b=-b ; end
    t0 = (1/a)*atan(b*tan((pi*a)/2));
    t=-t0:0.001:(pi/2); I=integrand2(t, x, a, t0); plot(t, I, 'r-')
end

function pdf = stab (x, a, b)
t0 = (1/a) * atan(b * tan((pi * a ) / 2));
if (x == 0) % Borak et al. (2005)
    pdf = gamma(1+(1/a)) * cos(t0) / pi / (1 + (b*tan((pi * a) / 2 ))^2)^(1/(2 * a));
else
    tol = 1e-9; display = 0;
    integ = quadv(@integrand , -t0, pi/2, tol, display, abs(x), a, t0);
    pdf = a * x^(1/(a-1)) / pi / abs(a-1) * integ;
end

function I = integrand (t, x, a, t0)
ct = cos(t); s = t0 + t;
v = (cos(a*t0))^(1/(a-1)) * (ct/sin(a*s))^(a / (a-1)) * (cos(a*s - t) / ct) ;
term = -x^(a/(a-1)); I = v.*exp(term * v);

function I = integrand2(t, x, a, t0)
ct = cos(t); s = t0 + t;
v = (cos(a*t0))^(1/(a-1)) * (ct / sin(a*s)).^(a / (a-1)).*(cos(a*s - t)./ ct);
term = -x.^(a/(a-1)); I = v.*exp(term * v);
