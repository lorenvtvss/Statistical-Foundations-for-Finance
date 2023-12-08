function I = fffplus(uvec, x, a, b, c, d, dopdf)
subs = 1; I = zeros(size(uvec));
for ii=1:length(uvec)
    u=uvec(ii);
    if subs==1, t=(1-u)/u; else t=u/(1-u); end
    if a==1, cf = exp(-c*abs(t)*(1+1i*b*(2/pi)*sign(t)*log(t)) + 1i*d*t);
    else
        cf=exp(-((c^a)*(abs(t))^a)*(1-1i*b*sign(t)*tan(pi*a/2)) + 1i*d*t);
    end
    z=exp(-1i*t*x).*cf; if dopdf==1,g=real(z); else g=imag(z)./t; end
    if subs==1, I(ii)=g*u^(-2); else I(ii)=g*(1-u )^(-2); end
end

