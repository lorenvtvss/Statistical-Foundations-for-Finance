function l = fff(uvec,x,a,b,dopdf)
subs=1; l=zeros(size(uvec));
for ii=1:length(uvec)
    u=uvec(ii);
    if subs==1, t=(1-u)/u; else t=u/(1-u); end
    if a==1; cf=exp(-abs(t)*(1+1i*b*(2/pi)*sign(t)*log(t)));
    else cf=exp(-((abs(t))^a)*(1-1i*b*sign(t)*tan(pi*a/2)));
    end
    z=exp(-1i*t*x).*cf;
    if dopdf==1, g=real(z); else g=imag(z)./t; end
    if subs==1, l(ii)=g*u^(-2); else l(ii)=g*(1-u)^(-2); end
end

