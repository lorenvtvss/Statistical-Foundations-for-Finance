function I=fff_generalized(uvec,x,a,b, c, d, dopdf, dogen)
    I = zeros(length(uvec), 1);
    for ii=1:length(uvec)
        u=uvec(ii);  t  = (1-u)/u;
        if dogen
            if a==1
                cf = exp(-c*abs(t)*(1 + 1i*b*(2/pi)*sign(t) * log(t)) + 1i*d*t);
            else
                cf = exp(-(c^a)*((abs(t))^a)*(1 - 1i*b*sign(t) * tan(pi*a/2))  + 1i*d*t);
            end
        else
            c=1; d=0;
            if a==1
                cf = exp( -abs(t)*( 1 + 1i*b*(2/pi)*sign(t) * log(t) ) ); 
            else
                cf = exp( - ((abs(t))^a) *( 1 - 1i*b*sign(t) * tan(pi*a/2) ) );
            end
        end
        z = exp(-1i*t*x) .* cf;
        if dopdf
            g=real(z);
        else
            g=imag(z)./t;
        end
        I(ii) = g / u^2;
    end
end