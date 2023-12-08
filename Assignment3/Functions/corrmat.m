function B = corrmat(d)
% function to create correlation matrix (positive definite)
r = -1 + (1+1)*rand(d);
M = tril(r,-1);
B = M + M' + eye(d);

while all(eig(B) > 0) == 0
    r = -1 + (1+1)*rand(d);
    M = tril(r,-1);
    B = M + M' + eye(d);
end
