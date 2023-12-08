function S = Stoy(cut, alpha, beta) 
if nargin<3, beta=0; end
cut= -cut; %weuseadifferentsignconvention
bbar = -sign(cut) * beta;
t0bar = (1 / alpha) * atan(bbar * tan(pi * alpha / 2));
% 'beta==0' => 'bbar==0' => 't0bar==0'
small = 1e-8; tol = 1e-8; abscut = abs(cut); display = 0; 
integ = quadl(@stoyint , -t0bar+small , pi/2-small , tol , ...
    display , abscut , alpha , t0bar ); 
S = alpha / (alpha-1) / pi * abscut * integ;

function I = stoyint(t, cut, a, t0bar)
s=t0bar + t;
g = sin(a*s - 2*t) ./ sin(a*s) - a * cos(t).^2 ./ sin(a*s).^2;
v = (cos(a*t0bar)).^(1 / (a - 1)) .* (cos(t) ./ sin(a*s)).^(a / (a-1)) ...
.* cos(a*s-t) ./ cos(t); term = -(abs(cut)^(a / (a-1))); I=g.*exp(term.*v);

