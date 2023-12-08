%%%% Assignment 1 %%%% 

%%% TODO
% Ex1: DONE
% Ex2: DONE
% Ex3: method4: conv --> done (optional: conv() with three variables) -->
% done; for-loop missing. 
% Ex4: DONE + optional: figure out how large the replication number (nb of
% samples) needs to be to ensure getting a certain number of correct
% significant digits, as compared to the theory values. 
% Ex5: DONE, just add the for loop for the different alphas 
% run all with higher sample size
% add descriptions and create latex document 

%% Exercise 1 - Plotting the pdf of stable: true density vs kernel density estimate. 


% question for Lore and Naina: do you prefer 

% We pick values of alpha=1.7, beta=-0.4, scale=2 and location=0.3. We also
% pick a reasonable sample size and we generate an array with first element
% as -20 and last element as 20, where each successive elements have a common 
% difference of 0.01. 
a=1.7;   % alpha
b=-0.4;  % beta 
c=2;     % scale
d=0.3;   % location 
n=1e6;  % sample size --> pick size first 1e4 and see how much it takes, until 1e6. 
xvec=-20:0.01:20;

% For creating the kernel desnsity estimate, we generate 1e6 simulated IID 
% variates from the S_{a,b}(d,c) distribution through the stabgen function.
% We afterwards simulate those IID stable variates through the ksdensity
% function, which returns a probability esimate f. 
% We finally plot the resulting density. 
eststab = stabgen(n,a,b,c,d,2);
[f,cd] = ksdensity(eststab, xvec);
figure, plot(cd, f, 'r-', 'linewidth', 2)

 
% The second line in our graph is the true density, based on analytic
% calculation through the asymstabplus function, which returns the pdf of
% the asymmetric stable. (calculating the actual theoretical values of a
% S_{a,b} distribution)
% We again plot the resulting density. 
truestab1 = asymstabplus(xvec,a,b,c,d);
hold on, plot(xvec, truestab1, 'b--', 'linewidth', 2)

% Adjusting the plot 
xlim([-20 20])
legend('Simulated PDF','Theoretical PDF', 'Location', 'NorthWest')
title('PDFs of the Stable Distribution')
xlabel("x"); ylabel("S_{1.7, -0.4}(0.3, 2)(x)")
set(gca, 'fontsize', 12)
hold off

saveas(gcf, 'Assignment1_ex1.png')

disp('Done')

% The goal was to plot the stable density using both the theoretical and
% simulation-based pdf (and a kernel density smoother). This is illustrated
% in Figure \ref{stablepdfcomp1}, with the red line showing the actual
% density, computed using numerical inversion of the characteristic
% function (codes shown in Program Listing ..?does he mean his code or 
% our functions?), and the blue line showing the kernel density estimate
% base on 1e6 (if we manage) replications (codes shown in listing ..(same
% question)). 
% For both we We pick values of alpha=1.7, beta=-0.4, scale=2 and location=0.3. We also
% pick a reasonable sample size and we generate an array with first element
% as -20 and last element as 20, where each successive elements have a common 
% difference of 0.01.


%% Exercise 2
a=1.7   %alpha
b1=-0.4 %beta 1
c1=2    %scale 1
d1=-0.5 %location 1
b2=-0.6 %beta 2
c2=4    %scale 2
d2=-0.3 %location 2
n = 1e6;
x=-30:0.01:30;

% since alpha is the same, it is easy to calculate beta, scale and location
% of the convolution
d=d1+d2;
c=((c1^a)+(c2^a))^(1/a);
b=((b1*c1^a)+(b2*c2^a))/((c1^a)+(c2^a));

% theoretical pdf of convolution --> ADJSUT ASYMSTAB!!
truestab_conv = asymstabplus(x,a,b,c,d);
figure, plot(x, truestab_conv, 'b-', 'linewidth', 2)

% kernel density estimate of simulated values of S
eststab = stabgen(n,a,b,c,d,2);
[f,cd] = ksdensity(eststab, x);
hold on, plot(cd, f, 'r--', 'linewidth', 2)

% adjusting the plot
title('Convolution of Stable Random Variables')
legend('True density','Kernel density')
xlim([-30 30])
xlabel("x"); ylabel("S_{a, b}(d, c)(x)");
set(gca, 'fontsize', 12)
hold off

saveas(gcf, 'Assignment1_ex2.png')

disp('Done')

%% Exercise 3 - Convolution of two independent stable random variables with different tail index alpha
a1=1.6;     %alpha 1   
%a1=1.5 
b1=0;               %beta 1
c1=1;               %scale 1
d1=0;               %location 1
a2=1.8;      %alpha 2    
%a2=1.9
b2=0;               %beta 2
c2=1;               %scale 2
d2=0;               %location 2

xvec1 = -15:.05:15;
xvec2 = -15:.05:15;


%%% Density convolution 1: using integration/convolution formula
conv1 = asymstabpdf_conv(xvec1, a1, b1, a2, b2, 1);


%%% Density convolution 2: using invsersion formula
conv2 = asymstab_invform(xvec1, a1, b1, a2, b2); 


%%% Density convolution 3: using simulation and kernel density
n = 1e6; 

% simulating random sample with different seeds to avoid correlation
random_conv_sim1 = stabgen(n,a1, b1, c1, d1, 5); 
random_conv_sim2 = stabgen(n,a2, b2, c2, d2, 265); 
random_conv_3 = random_conv_sim1 + random_conv_sim2; 
[f, c] = ksdensity(random_conv_3,xvec2); 


%%% Density convolution 4: using the conv() function
% Interval size for density convolution calculation through 
% the conv() function 
dx = 0.1; 
xvec= -15:dx:15;

% pdfs for the different alphas 
truestab_a1 = asymstabplus(xvec,a1,b1,c1,d1);
truestab_a2 = asymstabplus(xvec,a2,b2,c2,d2);

% Convolution calculation 
truestab_conv = conv(truestab_a1, truestab_a2, 'same')*dx; 


% plot
figure
plot(xvec1, conv1, 'b-', 'linewidth', 2)
hold on
plot(xvec1,conv2, 'r-.', 'linewidth', 2)
plot (c, f, 'g--', 'linewidth', 2)
plot (xvec',truestab_conv, 'm:', 'linewidth',2)
title('Comparison of Stable R.V. Convolution Methods')
legend({'Numeric Integration','Inversion Formula', 'Simulation', 'Matlab Conv Formula', '3 R.V. Conv Formula'}, 'FontSize',9)
xlim([-15 15])
xlabel("x"); ylabel("Sconv_{a, b}(d, c)(x)");
set(gca, 'fontsize', 12)
hold off

saveas(gcf, 'Assignment1_ex3_bonus.png')

disp('Done')

%% Exercise 3 - Bonus

% BONUS: density convolution using conv() function and three variables
a = 1.7;
b1 = 0; b2 = 0; b3 = 0; 
c1 = 1; c2 = 1; c3 = 1; 
d1 = 0; d2 = 0; d3 = 0; 
dx = 0.01; 
xvec= -15:dx:15;
n = 1e6;

% theoretical distribution
x = asymstabplus(xvec,a,b1,c1,d1);
y = asymstabplus(xvec,a,b2,c2,d2);
z = asymstabplus(xvec,a,b3,c3,d3); 

conv_bonus_1 = conv(x,y,'same')*dx
truestab_conv_bonus = conv(conv_bonus_1, z, 'same')*dx
figure, plot(xvec',truestab_conv_bonus, 'b-', 'linewidth',2)

% kernel density estimate of simulated values of S
b_conv = (b1 * c1^a + b2 * c2^a + b3 * c3^a)/(c1^a + c2^a + c3^a);
c_conv = (c1^a + c2^a + c3^a)^(1/a);
d_conv = d1 + d2 + d3;

eststab_conv3 = stabgen(n, a, b_conv, c_conv, d_conv, 8);
[f,cd] = ksdensity(eststab_conv3, xvec);
hold on, plot(cd, f, 'r--', 'linewidth', 2)

% adjusting plot
title('Convolution of three Stable R.V.')
legend('Simulated PDF','Theoretical PDF', 'Location', 'NorthWest')
xlim([-15 15])
xlabel("x"); ylabel("Sconv3_{a, b}(d, c)(x)");
set(gca, 'fontsize', 12)
hold off

saveas(gcf, 'Assignment1_ex3_bonus.png')


% apply the kernel density smoother

%% Exercise 4 
a = 1.7; 
b = 0; 
c = 0; 
d = 1; 
xi = 0.01; 
method = 1; 

% Theoretical ES based on a tail probability of xi=0.01 using the Stoyanov et al. result. 
[Theo_ES, Theo_VaR] = asymstableES(xi, a, b, c, d, method);
X = ['ES via Stoyanov et al: ', num2str(Theo_ES)]; 
disp(X);

% Simulation of the ES
nobs = 10^6;
P = stabgen(nobs, a, b, d, c, 0);
VaR = quantile(P, xi);
Plo = P(P <= VaR);
ES_simulated = mean(Plo); 
X = ['Simulated ES: ', num2str(ES_simulated)]; 
disp(X);

%% Exercise 5
a1=1.5;  %alpha 1    ->1.5
b1=0;    %beta 1
c1=1;    %scale 1
d1=0;    %location 1
a2=1.9;  %alpha 2    -> 1.9
b2=0;    %beta 2
c2=1;    %scale 2
d2=0;    %location 2

xi = [0.01 0.025 0.05];
xvec1 = -15:.05:15;

% simulated ES
n = 1e6;
X1 = stabgen(n, a1, b1, c1, d1, 1); X2 = stabgen(n, a2, b2, c2, d2, 2);
Sum = X1 + X2;
for xi=[0.01 0.025 0.05]
    VaR = quantile(Sum, xi);
    Plo = Sum(Sum <= VaR);
    ES_simulated = mean(Plo); 
    X = ['Simulated ES for xi=', num2str(xi), ': ', num2str(ES_simulated)]; 
    disp(X);
end

% Theoretical Stoyanov ES with estimated parameters
n = 1e6; method = 1;
X1 = stabgen(n, a1, b1, c1, d1, 1); X2 = stabgen(n, a2, b2, c2, d2, 2);
Sum = X1 + X2;
[alpha,beta,sigma,mu] = stablereg(Sum);
s1=['Alpha: ',num2str(alpha),' Beta: ',num2str(beta),' Sigma: ',num2str(sigma),' mu: ',num2str(mu)];
disp(s1)
for xi=[0.01 0.025 0.05]
    [Theo_ES, Theo_VaR] = asymstableES(xi, alpha, beta, mu, sigma, method);
    X2 = ['ES via Stoyanov et al. for xi=',num2str(xi),': ', num2str(Theo_ES)]; 
    disp(X2);
end

% We now repeat the above with aplha1=1.5 and alpha2=1.9


% -> the lower yi the more accurate?

% calc accuracy --> formula slides


