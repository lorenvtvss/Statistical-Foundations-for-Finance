function MLE = tlikmax0_modified(x,initvec,method)
tol=1e-5;
opts=optimset('Disp', 'none', 'LargeScale', 'Off', 'TolFun', tol, 'TolX', tol, 'Maxiter',200);
MLE = fminunc(@(param)tloglik_modified(param, x, method), initvec, opts);
end

function ll = tloglik_modified(param,x,method)
v=param(1);     % df
l=param(2);     % loc
c=param(3);     % scale
m=param(4);     % mu

% TBD: für jede method: ll berechnen (Aber was ist ll?)
if v<0.01, v=rand; end % An ad hoc way of prevent ing negative values
if c<0.01, c=rand; end % which works , but i s NOT recommended !

z = (x-l)/c;
if method == 'ddanct'
    %use stdnctpdfln_ j (is already log) (Try also with and without log)
    ll = stdnctpdfln_j(z, v, m)-log(c);
elseif method == 'matlab'
    %Duse matlab nct pdf fct & take log 
    ll = log(nctpdf(z, v, m))-log(c);
else disp('Invalid method used!')
end

ll = -sum(ll);
end
