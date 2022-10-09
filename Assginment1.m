%%%% Assignment 1 %%%%

%% Exercise 1 - Plotting density

a=1.7   % alpha
b=-0.4  % beta
c=2     % scale
d=0.3   % location
s=10^4  % sample size (10^6)

% 1. line: kernel density estimate based on, say, 1e6 (one million) simulated IID stable variates
% - simulate a stable random variable
x=-6:0.02:6;
x = stabgen(s,a,b,c,d,1)

% - pdf
f17 = asymstab(x,a,b)

% 2. line: true density, based on analytic calculation

% Final plot
figure
set(gca, 'fontsize', 16)
plot(x,f17,'r-','linewidth',2);
legend('\alpha = 1.7, \beta = -0.4','\alpha = ??, \beta = ??','Location','northwest')

%hold on
%plot(x,pdf2,'r-.');
%plot(x,pdf3,'k--');
%title('')
%xl = xlim;
%xlim([xl(1)-0.25 xl(2)])
%xlabel('')
%ylabel('')
%hold off


%% Exercise 2

%% Exercise 3

%% Exercise 4

%% Exercise 5
