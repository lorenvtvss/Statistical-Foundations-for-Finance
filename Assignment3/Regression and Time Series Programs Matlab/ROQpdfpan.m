function [f,F]=ROQpdfpan(rvec,A,B,N)
%[f,F]=ROQpdfpan(rvec,A,B,N=12)
%cdf and pdf of x'Ax/x'Bx, B>=0, x~N(0,I), at the values in rvec.
if nargin<4,N=12;end
T=length(A);F=rvec*0;f=F;
A=.5*(A+A');B=.5*(B+B');
yin=cos((2*[1:N]-1)*pi/(2*N));
for loop=1:length(rvec)
    r=rvec(loop);
    C=A-r*B;
    [P,L]=eig(C);
    wgt=diag(L);
    dwgt=-diag(P'*B*P);
    v=sum(wgt>0);
    sgn=(-1)^((v<=T/2)+1);
    [f(loop),F(loop)]=AS_R52(sgn*wgt,yin,dwgt);
    F(loop)=(sgn==-1)+sgn*F(loop);
end

function[f,F]=AS_R52(wgt,yin,dwgt)
[wgt,I]=sort(wgt,1,'descend');dwgt=dwgt(I);
m=length(wgt);v=sum(wgt>0);
F=1;f=0;
for j=1:ceil(0.5*v)
    Fj=mean(panintegrand(yin,wgt,j,dwgt,v),1);
    F=F+((-1)^j)*Fj(1);f=f+((-1)^j)*Fj(2);
end;

function[I]=panintegrand(tvec,wgt,j,dwgt,v)
m=length(wgt);
if j==(v+1)/2,
    outind=v;
    sv=1;
else
    outind=[2*j,2*j-1];
    sv=[1 ;-1];
end
inind=setdiff(1:m,outind);
dwgti=dwgt(inind);wgti=wgt(inind);
dwgto=dwgt(outind);wgto=wgt(outind);
I=tvec*0;
for tloop=1:length(tvec)
    t=tvec(tloop);
    yv=sum(wgto)-wgto'*sv*t;
    I(tloop,1)=exp(-.5*sum(log(abs(1-2*wgti/yv))));
    dyv=sum(dwgto)-dwgto'*sv*t;
    gprime=sum((dwgti-dyv/yv*wgti)./(yv-2*wgt(inind)));
    I(tloop,2)=I(tloop)*gprime;
end;
