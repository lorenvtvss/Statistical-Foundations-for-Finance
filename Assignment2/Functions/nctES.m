function [ES,VaR] = nctES(xi,v,theta)
howfar = nctinv(1e-8,v,theta); %howfarintothelefttailtointegrate
VaR = nctinv(xi,v,theta); %matlabroutineforthequantile
I = quadl(@int,howfar,VaR,1e-6,[],v,theta);ES=I/xi;
end

function I=int(u,v,theta);
pdf=nctpdf(u,v,theta);
I=u.*pdf;
end
