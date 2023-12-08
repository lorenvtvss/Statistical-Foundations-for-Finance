function [f,F] = asymstab_generalized(xvec, a, b, c, d)
    bordertol=1e-8; lo=bordertol; hi=1-bordertol; tol=1e-6;
    xl=length(xvec); F=zeros(xl,1); f=F;
    if nargin>3
        dogen=1;
        for loop=1:length(xvec)
            x=xvec(loop); % dopdf=1;
            f(loop) = quadl(@fff_generalized, lo, hi, tol, [], x, a, b, c, d, 1, dogen)/pi;
            if nargout>1
                F(loop) = 0.5 - (1/pi)*quadl(@fff_generalized, lo, hi, tol, [], x, a, b, c, d, 0, dogen);
            end
        end
    else
        dogen=0;
        for loop=1:length(xvec)
            x=xvec(loop); % dopdf=1;
            f(loop)= quadl(@fff_generalized,lo,hi,tol,[],x,a,b,1,0,1,dogen) / pi;
            if nargout>1
                F(loop)=0.5-(1/pi)* quadl(@fff_generalized,lo,hi,tol,[],x,a,b,1,0,0,dogen);
            end
        end
    end
end