function MLE = tlikmax0(x,initvec)
tol=1e-5;
opts=optimset('Disp', 'none', 'LargeScale', 'Off', 'TolFun', tol, 'TolX', tol, 'Maxiter',200);
MLE = fminunc(@(param)tloglik(param, x), initvec, opts);
end

function ll = tloglik(param,x)
v=param(1); mu=param(2); c=param(3);
if v<0.01, v=rand; end % An ad hoc way of prevent ing negative values
if c<0.01, c=rand; end % which works , but i s NOT recommended !
K=beta(v/2,0.5)*sqrt(v); z=(x-mu)/c;
ll = -log(c)-log(K)-((v+1)/2)*log(1+(z.^2)/v); 
ll = -sum(ll);
end
