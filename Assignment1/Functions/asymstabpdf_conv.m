function f = asymstabpdf_conv(xvec1, a1, b1, a2, b2, plotintegrand)
% pdf of convolution of two asymmetric stable with different alphas

if nargin<4 , plotintegrand =0; end
xl=length(xvec1); f=zeros(xl, 1); n = length(xvec1);

for loop=1:n
    x=xvec1(loop);
    f(loop)=integral(@(xvec2) integrand_conv(xvec2,x,a1,b1,a2,b2), -Inf, Inf);
end





