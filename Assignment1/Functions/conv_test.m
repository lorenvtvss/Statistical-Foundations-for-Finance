%% convolution example
alpha1 = 1.6;
alpha2 = 1.8;
beta = 0;
c = 1;
d = 0;

% size of intervals to calculate the density using conv function
dx = 0.01;
% x to evaluate
x = -10:dx:10;

% calculate range for the convolution
%xconv = linspace(2*x(1), 2*x(end), 2*length(x)-1);

% generate pdf's of each Stable
s1 = asymstab_generalized(x, alpha1, beta, c, d);
s2 = asymstab_generalized(x, alpha2, beta, c, d);

% calculate the convolution
sconv = conv(s1, s2, 'same')*dx;

% plot
plot(x', sconv, 'b', 'LineWidth', 2)