N = 3; 
d = 1; % The diagonal values
t = triu(bsxfun(@min,d,-d.').*rand(N),1); % The upper trianglar random values
M = diag(d)+t+t.'; % Put them together in a symmetric matrix
disp(M)
disp(eig(M)); 

nu = 4; 
mu = 0; 
scale = 1; 
n = 3; 
m = 3; 
x1 = ((-3:.2:3) - mu)/scale; x2 = ((-3:.2:3) - mu)/scale; x3 = eye(n,m); 
[X1,X2,X3] = meshgrid(x1,x2,x3);
F = mvtpdf([X1(:) X2(:) X3(:)],M,nu);
disp(F)
F = reshape(F,[],length(x2),length(x1));
surf(x1,x2,z);
caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
axis([-3 3 -3 3 0 .2])
xlabel('x1'); ylabel('x2'); zlabel('Probability Density');