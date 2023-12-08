function f = asymstab_invform(xvec, a1, b1, a2, b2)

bordertol=1e-8; lo=bordertol; hi=1-bordertol; tol=1e-7;
xl=length(xvec); F=zeros(xl,1); f=F;

for loop=1:length(xvec)
    x=xvec(loop); 
    f(loop)=integral(@(uvec) fff_inv(uvec,x,a1,b1,a2,b2,1), lo, hi) / pi;
end