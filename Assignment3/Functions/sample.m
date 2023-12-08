% This function generates samples of a multivariate Student-t distribution.
function X=sample(n,mu,nu,sigma)
[a,b]=eig(sigma);
sqrt_sigma=a'*(b.^.5)*a;
d=size(mu,1);
X=repmat(mu,1,n)+(sqrt_sigma*normrnd(0,1,[d,n]))./repmat((gamrnd(nu/2,2/nu,[1,n])).^.5,d,1);
end