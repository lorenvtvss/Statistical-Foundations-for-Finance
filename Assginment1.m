%%%% Assignment 1 %%%%

%% Exercise 1 - Plotting üdf of stable

a=1.7;   % alpha
b=-0.4;  % beta 
c=2;     % scale
d=0.3;   % location
n=4000;  % sample size
x=-20:.001:20;  %x=-6:0.02:6;

%%% kernel density estimate 
% generating a random sample of size n from the S_{a,b}(d, c) distribution
% and plot the resulting density
eststab = stabgen(n,a,b,c,d,2)
[f,c] = ksdensity(eststab, x);
figure, plot(c, f, 'r--', 'linewidth', 2)
xlim([-20 20])


%%% true density 
% calculating the actual theoretical values of a S_{a,b} distribution
truestab = asymstab(x, a, b);
%truestab = asymstabplus(x, a, b, c, d);
hold on, plot(x, truestab, 'b-', 'linewidth', 2), hold off

% adjust plot
legend('Simulated PDF', ...
       'Theoretical PDF', 'Location', 'NorthWest')
title('PDFs of Stable Distribution')
xlabel("x"); ylabel("S_{1.7, -0.4}(2, 0.3)(x)")
set(gca, 'fontsize', 16)


%% Exercise 2

%% Exercise 3

%% Exercise 4

%% Exercise 5
