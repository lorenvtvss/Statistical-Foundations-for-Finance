function l = fff_inv(uvec,x,a1,b1,a2,b2,dopdf)

subs=1; l = zeros(size(uvec));
for ii=1:length(uvec)
    u=uvec(ii);
    if subs==1, t =(1-u)/u ; else t=u/(1-u) ; end
    cf1 = exp(-((abs(t))^a1)*(1-1i*b1*sign(t)*tan(pi*a1/2)));
    cf2 = exp(-((abs(t))^a2)*(1-1i*b2*sign(t)*tan(pi*a2/2)));
    cf = cf1.*cf2;
    z=exp(-1i*t*x).*cf;
    if dopdf==1, g=real(z); else g=imag(z)./t; 
    end
    l(ii)=g*u^(-2); 
end