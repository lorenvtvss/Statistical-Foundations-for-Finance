function [param,stderr,iters,loglik,Varcov] = tlikmax(x,initvec)
%%%%%%%%    df      mu  c
bound.lo=   [1      1   0.01];
bound.hi=   [100    1   100 ];
bound.which=[1      0   1   ];
%Inthiscase,asbound.whichformuiszero,muwillnotbe
%restricted.Assuch,thevaluesfor.loand.hiareirrelevant
maxiter=100; tol=1e-3;%changetheseasyouseefit
opts=optimset('Display','notify-detailed','Maxiter',maxiter,'TolFun',tol,'TolX',tol,'LargeScale','Off');
[pout,fval,exitflag,theoutput,grad,hess]=fminunc(@(param)tloglik(param,x,bound),einschrk(initvec,bound),opts);
V=inv(hess); %Don'tnegate:weworkwiththenegoftheloglik
[param,V]=einschrk(pout,bound,V); %Transformback,applydeltamethod
param=param'; Varcov=V;
stderr=sqrt(diag(V)); %Approxstderroftheparams
loglik=-fval; %Thevalueoftheloglikatitsmaximum.
iters=theoutput.iterations; %Numberofloglikfunctionevals
end


function ll = tloglik(param,x,bound)
if nargin<3,bound=0; end
if isstruct(bound),paramvec=einschrk(real(param),bound,999);
else paramvec=param; end

v=paramvec(1); mu=paramvec(2); c=paramvec(3);
K=beta(v/2,0.5)*sqrt(v); z=(x-mu)/c;
ll=-sum(-log(c)-log(K)-((v+1)/2)*log(1+(z.^2)/v));
end

function [pout,Vout] = einschrk(pin,bound,Vin)
lo=bound.lo; hi=bound.hi ; welche=bound.which;
if nargin<3
    trans=sqrt((hi-pin)./(pin-lo));pout=(1-welche).*pin+welche.*trans;
    Vout=[];
else
    trans=(hi+lo.*pin.^2)./(1+pin.^2); pout=(1-welche).*pin+welche.*trans;
    %nowadjustthestandarderrors
    trans=2*pin.*(lo-hi)./(1+pin.^2).^2;
    d=(1-welche)+welche.*trans; %eitherunityordeltamethod.
    J=diag(d); Vout=J*Vin*J;
end
end
