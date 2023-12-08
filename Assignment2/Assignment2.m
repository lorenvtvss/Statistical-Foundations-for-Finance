use%%%% Assignment 2 %%%%%


%% Question 1:
%%%% build upon student_t.m script

% Fixing the 3 true parameters: tailprob alpha, dof, left tail quantile of location 0, scale 1 distribution
alpha = 0.1; df = 4; c01 = tinv(alpha , df); % alpha = 0.01

% Parameters 
loc = 0; scale = 1;         % loc = 1; scale = 2;
rep = 100;                  % 1000 Repititions
T_samples = [50, 100];          % [250 500 2000] Sample size
B = 100;                    % 1000 Bootstrap samples
rand('twister',6)           % set seed value so we can replicate results
initvec = [df loc scale];   % for parametric MLE

% Initializing vectors
ESvec_nonparam = zeros(B, 1); 
ESvec_param_book = zeros(B, 1); 
ESvec_param_matlab = zeros(B, 1); 
ci_length_nonparam = zeros(rep, length(T_samples)); % length of ci interval for each rep
ci_length_param_book = zeros(rep, length(T_samples)); 
ci_length_param_matlab = zeros(rep, length(T_samples)); 
ci_trueES_nonparam = zeros(rep, length(T_samples)); % whether interval contains trueES for each rep
ci_trueES_param_book = zeros(rep, length(T_samples)); 
ci_trueES_param_matlab = zeros(rep, length(T_samples)); 


% True ES for student t
truec = loc+scale*c01;                                          % left tail quantile c
ES_01_analytic = -tpdf(c01,df)/tcdf(c01,df)*(df+c01^2)/(df-1);
trueES = loc+scale*ES_01_analytic;                              % true theoretical ES


% Simulate "rep" repetitions of an IID T-length sequence of Student t
% For each rep claculate bootstrap CI 90% based on B bootstrap replications
for t = 1:length(T_samples)
    for r = 1:rep
        %%% generating data points from student t distribution
        data=loc+scale*trnd(df,T_samples(t),1);               
        
        
        %%% parametric bootstrap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BOOK
        % book code -> Listing 2.9, p.56 
        %y = trnd(df,T_samples(t),B)+1;
        %phatvec = 1./mean(trnd(df,T_samples(t),B)+1) ; % MLE
        %ci_param_book = quantile(phatvec, [alpha/2 1-alpha/2]); % parametric
        
        % book code -> Listing 4.5, p.139
        mle_param_book = tlikmax0(data, initvec);

        % Now that I have the parameters' estimation through MLE based on the bootstrap data set
        % I want to compute the theoretical ES based on those parameter values. 
        for b=1:B
            param_bs_sample_book = mle_param_book(2)+mle_param_book(3)*trnd(mle_param_book(1),T_samples(t),1);
            VaR_param_book = quantile(param_bs_sample_book, alpha); 
            temp_book = param_bs_sample_book(param_bs_sample_book<=VaR_param_book); 
            ESvec_param_book(b) = mean(temp_book);
        end 

        % length of CI (parametric bs)
        ci_param_book = quantile(ESvec_param_book,[alpha/2 1-alpha/2]); 
        %ci_param_book = quantile(phatvec, [alpha/2 1-alpha/2]);
        low_param_book = ci_param_book(1); high_param_book = ci_param_book(2); 
        ci_length_param_book(r,t) = high_param_book-low_param_book; 

        % whether or not the interval contains the TRUE ES (parametric bs)
        ci_trueES_param_book(r,t) = (trueES>low_param_book)&(trueES<high_param_book); 

        
        
        % MATLAB BUILT IN:
        mle_param_matlab = mle(data, 'Distribution', 'tLocationScale'); % [loc,scale,nu] output
        for b=1:B
            param_bs_sample_matlab = mle_param_matlab(1)+mle_param_matlab(2)*trnd(mle_param_matlab(3),T_samples(t),1); 
            VaR_param_matlab = quantile(param_bs_sample_matlab, alpha); 
            temp = param_bs_sample_matlab(param_bs_sample_matlab<=VaR_param_matlab);
            ESvec_param_matlab(b) = mean(temp);
        end
        
        % length of CI (parametric bs)
        ci_param_matlab = quantile(ESvec_param_matlab, [alpha/2 1-alpha/2]);
        low_param_matlab = ci_param_matlab(1); high_param_matlab = ci_param_matlab(2);
        ci_length_param_matlab(r,t) = high_param_matlab-low_param_matlab;
        
        % whether or not the interval contains the TRUE ES (parametric bs)
        ci_trueES_param_matlab(r,t) = (trueES>low_param_matlab)&(trueES<high_param_matlab);
        %disp(ci_trueES_nonparam(r,t));
        
        
        
        %%% non-parametric bootstrap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % book code -> Listing 2.9, p.56
        for b=1:B
            ind = unidrnd(T_samples(t),[T_samples(t),1]);
            nonparam_bs_sample=data(ind);
            %phatvec_nonpara_bs(b) = 1/mean(nonparam_bs); % necessary?
            VaR_nonpara = quantile(nonparam_bs_sample, alpha); 
            temp_nonparam = nonparam_bs_sample(nonparam_bs_sample<=VaR_nonpara);
            ESvec_nonparam(b) = mean(temp_nonparam);
        end

        %%% length of CI (non-parametric bs)
        ci_nonparam = quantile(ESvec_nonparam, [alpha/2 1-alpha/2]);
        low_nonparam = ci_nonparam(1); high_nonparam = ci_nonparam(2);
        ci_length_nonparam(r,t) = high_nonparam-low_nonparam;
        
        %%% whether or not the interval contains the TRUE ES (non-parametric bs)
        ci_trueES_nonparam(r,t) = (trueES>low_nonparam)&(trueES<high_nonparam);
        %disp(ci_trueES_nonparam(r,t));
        
        
    end
    
    
    % Printing some stuff
    disp('================================')
    X = ['True ES:', num2str(trueES)]; 
    disp(X);

    text = ['For T = ', int2str(T_samples(t)),' Rep = ',int2str(rep),' B = ',int2str(B)];
    disp(text)
    disp('--------------------------------')
    disp('Parametric (Book codes):')
    mean_ci_length_param_book = mean(ci_length_param_book(:,t));
    A1 = ['Mean CI length:', num2str(mean_ci_length_param_book)];
    disp(A1)
    sum_bool1 = sum(ci_trueES_param_book(:,t));
    coverage_ratio1 = sum_bool1/length(ci_trueES_param_book(:,t));
    B1 = ['Coverage ratio:', num2str(coverage_ratio1)]; 
    disp(B1)

    disp('--------------------------------')
    disp('Parametric (MATLAB built-in):')
    mean_ci_length_param_matlab = mean(ci_length_param_matlab(:,t));
    A2 = ['Mean CI length:', num2str(mean_ci_length_param_matlab)];
    disp(A2)
    sum_bool2 = sum(ci_trueES_param_matlab(:,t));
    coverage_ratio2 = sum_bool2/length(ci_trueES_param_matlab(:,t));
    B2 = ['Coverage ratio:', num2str(coverage_ratio2)]; 
    disp(B2)

    disp('--------------------------------')
    disp('Non-parametric:')
    mean_ci_length_nonparam = mean(ci_length_nonparam(:,t));
    A3 = ['Mean CI length:', num2str(mean_ci_length_nonparam)]; 
    disp(A3)
    sum_bool3 = sum(ci_trueES_nonparam(:,t));
    coverage_ratio3 = sum_bool3/length(ci_trueES_nonparam(:,t));
    B3 = ['Coverage ratio:', num2str(coverage_ratio3)]; 
    disp(B3)
    
    % additional check on the parameters 
    disp('--------------------------------')
    Z = ['BS parameters book MLE:', num2str(mle_param_book)]; 
    disp(Z)
    Y = ['BS parameters matlab MLE:', num2str(mle_param_matlab)]; 
    disp(Y)




%     % Now make a nice, visually appealing graphic:
%     actualcoverage=mean(ci_trueES);
%     disp(actualcoverage);
%     figure
%     histogram(ESvec), ax=axis;
%     set(gca,'fontsize',12)
%     line ([ trueES trueES ] ,[0 ax(4)], 'color', 'g ', 'linewidth',3)
%     xlabel('ES value (simulation and true as vertical line)')
%     title(['Simulated Stud t Empirical ES, T=',int2str(T_samples(t)),' obs'])
end






%% Question 2: 

% Same procedure applied in Q1, but such that the true data are
% non-central Student t, but wrongly assuming a 'regular' location-scale
% Student t. 
% This means, your codes for the parametric and nonparametric bootstrap 
% stay the same, but you USE AS THE TRUE DATA simulations from the 
% noncentral student t.

%%% PART 1: simulation of a noncentral t (NCT)

% Book code and function --> does not work
mu = [0 -1 -2 -3];  
df = 6; % df = 3         
x = (-15:0.1:15)'; 

for j=1:length(mu)
    %NCT_sim = stdnctpdfln_j(x, df, mu); 
    NCT_sim = exp(stdnctpdfln_j(x, df, mu(j))); 

    % MATLAB built in function --> works, but we have to compare both 
    %nct = log(nctpdf(x, df, mu)); 
    nct = nctpdf(x, df, mu(j)); 

    hold on, 
    plot(x, NCT_sim, 'b','LineWidth', 2)
    plot(x, nct, 'r--','LineWidth', 2)
    xlim([-15 15])
    legend('d.d.a. NCT','MATLAB built-in fct.')
    title('Simulated PDF of the NCT')
    xlabel("x"); ylabel("NCT(\nu,\gamma)(x)")
    set(gca, 'fontsize', 12)

    name = ['Assignment2_ex2_df',int2str(df),'.png'];

    saveas(gcf,name)
    hold off
end



%%% PART 2: compute the true ES of NCT distribution - first via simulation,
% and then by using the integral definition of the NCT (namely solve for
% the alpha quantile, and then do numeric integration of the NCT density 
% multiplied by x). 
clear
clc
mu = -1;
df = 6;
x = (-15:0.1:15)'; 

% ES via simulation
alpha = 0.1; loc=0; scale=1; n = 10^6;
%nct_rnd_sample = loc+scale*nctrnd(df,mu,100,1);
norm=normrnd(mu,1,n,1); chi2=chi2rnd(df,n,1);         % central chi2
%norm=normrnd(mu,1, n, 1); chi2=ncx2rnd(df, 0, n, 1);  % noncentral chi2
nct_rnd_sample = loc+scale*(norm./sqrt(chi2/df));
VaR = quantile(nct_rnd_sample, alpha);
temp_2 = nct_rnd_sample(nct_rnd_sample <= VaR);
ES_simulated = mean(temp_2); 
disp(['Simulated ES: ', num2str(ES_simulated)]);

% ES using integral definition of NCT (Ila)
c01_nct = nctinv(alpha , df, mu); 
ES_01_nct = @(x) x.*nctpdf(x, df, mu); 
ES_nct_int= integral(ES_01_nct, -Inf, c01_nct)/alpha; 
disp(['ES via Numeric Integration:', num2str(ES_nct_int)]);



%%% PART 3: your codes for the parametric and nonparametric bootstrap 
% stay the same, but you USE AS THE TRUE DATA simulations from the 
% noncentral student t.


% Fixing the 3 true parameters: tailprob alpha, dof, left tail quantile of location 0, scale 1 distribution
%alpha = 0.1; df = 4; c01 = tinv(alpha , df); % alpha = 0.01

% Parameters 
loc = 0; scale = 1;         % loc = 1; scale = 2;
rep = 1000;                  % 1000 Repititions
T_samples = [500];           % [100 500 2000] Sample size
B = 1000;                    % 1000 Bootstrap samples
df = df;                     % df = 3, 6
rand('twister',6)           % set seed value so we can replicate results
initvec = [df loc scale];   % for parametric MLE

% Initializing vectors
ESvec_nonparam = zeros(B, 1); 
ESvec_param_book = zeros(B, 1); 
ESvec_param_matlab = zeros(B, 1); 
ci_length_nonparam = zeros(rep, length(T_samples)); % length of ci interval for each rep
ci_length_param_book = zeros(rep, length(T_samples)); 
ci_length_param_matlab = zeros(rep, length(T_samples)); 
ci_trueES_nonparam = zeros(rep, length(T_samples)); % whether interval contains trueES for each rep
ci_trueES_param_book = zeros(rep, length(T_samples)); 
ci_trueES_param_matlab = zeros(rep, length(T_samples)); 


% True ES for student t
trueES = ES_nct_int;
%trueES = ES_simulated;


% Simulate "rep" repetitions of an IID T-length sequence of Student t
% For each rep claculate bootstrap CI 90% based on B bootstrap replications
for t = 1:length(T_samples)
    for r = 1:rep
        %%% generating data points from student t distribution
        norm=normrnd(mu,1, T_samples(t), 1); chi2=chi2rnd(df,T_samples(t),1);  % noncentral chi2
        data = loc+scale*(norm./sqrt(chi2/df));
        
        
        %%% parametric bootstrap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BOOK
        % book code -> Listing 2.9, p.56 
        %y = trnd(df,T_samples(t),B)+1;
        %phatvec = 0/mean(trnd(df,T_samples(t),B)+1) ; % MLE
        %ci_param_book = quantile(phatvec, [alpha/2 1-alpha/2]); % parametric
        
        % book code -> Listing 4.5, p.139
        mle_param_book = tlikmax0(data, initvec);

        % Now that I have the parameters' estimation through MLE based on the bootstrap data set
        % I want to compute the theoretical ES based on those parameter values. 
        for b=1:B
            param_bs_sample_book = mle_param_book(2)+mle_param_book(3)*trnd(mle_param_book(1),T_samples(t),1);
            VaR_param_book = quantile(param_bs_sample_book, alpha); 
            temp_book = param_bs_sample_book(param_bs_sample_book<=VaR_param_book); 
            ESvec_param_book(b) = mean(temp_book);
        end 

        % length of CI (parametric bs)
        ci_param_book = quantile(ESvec_param_book,[alpha/2 1-alpha/2]); 
        %ci_param_book = quantile(phatvec, [alpha/2 1-alpha/2]);
        low_param_book = ci_param_book(1); high_param_book = ci_param_book(2); 
        ci_length_param_book(r,t) = high_param_book-low_param_book; 

        % whether or not the interval contains the TRUE ES (parametric bs)
        ci_trueES_param_book(r,t) = (trueES>low_param_book)&(trueES<high_param_book); 

        
        
        % MATLAB BUILT IN:
        mle_param_matlab = mle(data, 'Distribution', 'tLocationScale'); % [loc,scale,nu] output
        for b=1:B
            param_bs_sample_matlab = mle_param_matlab(1)+mle_param_matlab(2)*trnd(mle_param_matlab(3),T_samples(t),1); 
            VaR_param_matlab = quantile(param_bs_sample_matlab, alpha); 
            temp = param_bs_sample_matlab(param_bs_sample_matlab<=VaR_param_matlab);
            ESvec_param_matlab(b) = mean(temp);
        end
        
        % length of CI (parametric bs)
        ci_param_matlab = quantile(ESvec_param_matlab, [alpha/2 1-alpha/2]);
        low_param_matlab = ci_param_matlab(1); high_param_matlab = ci_param_matlab(2);
        ci_length_param_matlab(r,t) = high_param_matlab-low_param_matlab;
        
        % whether or not the interval contains the TRUE ES (parametric bs)
        ci_trueES_param_matlab(r,t) = (trueES>low_param_matlab)&(trueES<high_param_matlab);
        %disp(ci_trueES_nonparam(r,t));
        
        
        
        %%% non-parametric bootstrap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % book code -> Listing 2.9, p.56
        for b=1:B
            ind = unidrnd(T_samples(t),[T_samples(t),1]);
            nonparam_bs_sample=data(ind);
            %phatvec_nonpara_bs(b) = 1/mean(nonparam_bs); % necessary?
            VaR_nonpara = quantile(nonparam_bs_sample, alpha); 
            temp_nonparam = nonparam_bs_sample(nonparam_bs_sample<=VaR_nonpara);
            ESvec_nonparam(b) = mean(temp_nonparam);
        end

        %%% length of CI (non-parametric bs)
        ci_nonparam = quantile(ESvec_nonparam, [alpha/2 1-alpha/2]);
        low_nonparam = ci_nonparam(1); high_nonparam = ci_nonparam(2);
        ci_length_nonparam(r,t) = high_nonparam-low_nonparam;
        
        %%% whether or not the interval contains the TRUE ES (non-parametric bs)
        ci_trueES_nonparam(r,t) = (trueES>low_nonparam)&(trueES<high_nonparam);
        %disp(ci_trueES_nonparam(r,t));
        
        
    end
    
    
    % Printing some stuff
    disp('================================')
    X = ['True ES:', num2str(trueES)]; 
    disp(X);

    text = ['For T = ', int2str(T_samples(t)),' Rep = ',int2str(rep),' B = ',int2str(B)];
    disp(text)
    disp('--------------------------------')
    disp('Parametric (Book codes):')
    mean_ci_length_param_book = mean(ci_length_param_book(:,t));
    A1 = ['Mean CI length:', num2str(mean_ci_length_param_book)];
    disp(A1)
    sum_bool1 = sum(ci_trueES_param_book(:,t));
    coverage_ratio1 = sum_bool1/length(ci_trueES_param_book(:,t));
    B1 = ['Coverage ratio:', num2str(coverage_ratio1)]; 
    disp(B1)

    disp('--------------------------------')
    disp('Parametric (MATLAB built-in):')
    mean_ci_length_param_matlab = mean(ci_length_param_matlab(:,t));
    A2 = ['Mean CI length:', num2str(mean_ci_length_param_matlab)];
    disp(A2)
    sum_bool2 = sum(ci_trueES_param_matlab(:,t));
    coverage_ratio2 = sum_bool2/length(ci_trueES_param_matlab(:,t));
    B2 = ['Coverage ratio:', num2str(coverage_ratio2)]; 
    disp(B2)

    disp('--------------------------------')
    disp('Non-parametric:')
    mean_ci_length_nonparam = mean(ci_length_nonparam(:,t));
    A3 = ['Mean CI length:', num2str(mean_ci_length_nonparam)]; 
    disp(A3)
    sum_bool3 = sum(ci_trueES_nonparam(:,t));
    coverage_ratio3 = sum_bool3/length(ci_trueES_nonparam(:,t));
    B3 = ['Coverage ratio:', num2str(coverage_ratio3)]; 
    disp(B3)


% additional check on the parameters 
    disp('--------------------------------')
    Z = ['BS parameters book MLE:', num2str(mle_param_book)]; 
    disp(Z)
    Y = ['BS parameters matlab MLE:', num2str(mle_param_matlab)]; 
    disp(Y)

end




%% Question 3 

% Repeat Q2 but using SYMMETRIC stable Paretian distribution. 
% Possible to use MATLAB's built in function for it. 

% Parameters
a1 = 1.6; %a1 = 1.8; 
b = 0; c = 1; d = 0; n = 1e7; 
xvec = -20:0.01:20; 
alpha = 0.1; 
method = 1; 
seed = 2;
%rng(seed, 'Twister')

T = 200;

loc = 0; scale = 1;         
rep = 1000;                  % 1000 Repititions
alpha = 0.1; 
df = 6; % df = 6;
T_samples = 100;            % [100 500 2000] Sample size
B = 1000;                    % 1000 Bootstrap samples
                   % set seed value so we can replicate results
%initvec = [df loc scale];   % for parametric MLE
initvec = [df loc scale];

%%% PART 1: simulation of a symmetric stable (meaning beta = 0) Paretian 
% distribution using alpha1 = 1.6 and alpha2 = 1.8. 


% Creating the kernel density estimate by generating 1e6 simulated IID
% variates from the stable distribution and simulating them. 
eststab = stabgen(n,a1,b,c,d, seed); 
%eststab = stabgen(n,a1,b,c,d);
[f,cd] = ksdensity(eststab,xvec); 
figure, hold on, plot(cd, f , 'r-', 'linewidth', 2)

t_dist = trnd(df,T,1);
[f2,cd2] = ksdensity(t_dist,xvec); 
plot(cd2, f2 , 'b--', 'linewidth', 2)

% adjusting the plot
xlim([-20 20])
legend('Simulated PDF of the Symmtric Stable Paretian distribution', 'Location', 'NorthWest')
title('Simulated PDF of the Symmetric Stable Paretian distribution')
hold off

%%% PART 2: compute the true ES of the symmetric stable distribution - 
% both via the Stoyanov et al. result and by simulation. 

% Theoretical ES based on a tail probability xi=0.01 using the Stoyanov et
% al. result
[Theo_ES, Theo_VaR] = asymstableES(alpha , a1, b, d, c, method); 
X = ['ES via Stoyanov et al: ', num2str(Theo_ES)];
disp(X);

% Simulation of ES
nobs = 10^6;
P = stabgen(nobs, a1, b, c, d, seed); 
VaR = quantile(P, alpha);
Plo=P(P<=VaR);
ES_simulated = mean(Plo);
X = ['Simulated ES: ' , num2str(ES_simulated)]; 
disp(X);

%%% PART 3: your codes for the parametric and nonparametric bootstrap 
% stay the same, but you USE AS THE TRUE DATA simulations from the 
% symmetric Stable Paretian distribution.


% Fixing the 3 true parameters: tailprob alpha, dof, left tail quantile of location 0, scale 1 distribution
%alpha = 0.1; df = 4; c01 = tinv(alpha , df); % alpha = 0.01


% Initializing vectors
ESvec_nonparam = zeros(B, 1); 
ESvec_param_book = zeros(B, 1);
ESvec_param_matlab = zeros(B, 1); 
ci_length_nonparam = zeros(rep, length(T_samples)); % length of ci interval for each rep
ci_length_param_book = zeros(rep, length(T_samples)); 
ci_length_param_matlab = zeros(rep, length(T_samples)); 
ci_trueES_nonparam = zeros(rep, length(T_samples)); % whether interval contains trueES for each rep
ci_trueES_param_book = zeros(rep, length(T_samples)); 
ci_trueES_param_matlab = zeros(rep, length(T_samples)); 


% True ES for student t
trueES = Theo_ES;
%trueES = ES_simulated;


% Simulate "rep" repetitions of an IID T-length sequence of Student t
% For each rep claculate bootstrap CI 90% based on B bootstrap replications
for t = 1:length(T_samples)
    for r = 1:rep
        %%% generating data points from student t distribution
        % norm=normrnd(mu,1, T_samples(t), 1); chi2=ncx2rnd(df, 0, T_samples(t), 1);  % noncentral chi2
        % data = loc+scale*(stabgen(n,a1,b,c,d,2));
        % data = loc+scale*(ksdensity(eststab,xvec));
        % data = c+d*(eststab);
        % data = loc+scale*(stabgen(nobs,a1,b,c,d,seed)); 
        % data = stabgen(T_samples(t),a1,b,d,c,seed); 
        % data = ksdensity(eststab,xvec); 
        data = stabgen(T_samples(t), a1, 0, scale, loc); 

        %%% parametric bootstrap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BOOK
        % book code -> Listing 2.9, p.56 
        %y = trnd(df,T_samples(t),B)+1;
        %phatvec = 1./mean(trnd(df,T_samples(t),B)+1) ; % MLE
        %ci_param_book = quantile(phatvec, [alpha/2 1-alpha/2]); % parametric
        
        % book code -> Listing 4.5, p.139
        mle_param_book = tlikmax0(data, initvec);

        % Now that I have the parameters' estimation through MLE based on the bootstrap data set
        % I want to compute the theoretical ES based on those parameter values. 
        for b=1:B
            param_bs_sample_book = mle_param_book(2)+mle_param_book(3)*trnd(mle_param_book(1),T_samples(t),1);
            VaR_param_book = quantile(param_bs_sample_book, alpha); 
            temp_book = param_bs_sample_book(param_bs_sample_book<=VaR_param_book); 
            ESvec_param_book(b) = mean(temp_book);
        end 

        % length of CI (parametric bs)
        ci_param_book = quantile(ESvec_param_book,[alpha/2 1-alpha/2]); 
        %ci_param_book = quantile(phatvec, [alpha/2 1-alpha/2]);
        low_param_book = ci_param_book(1); high_param_book = ci_param_book(2); 
        ci_length_param_book(r,t) = high_param_book-low_param_book; 

        % whether or not the interval contains the TRUE ES (parametric bs)
        ci_trueES_param_book(r,t) = (trueES>low_param_book)&(trueES<high_param_book);
        
        % MATLAB BUILT IN:
        mle_param_matlab = mle(data, 'Distribution', 'tLocationScale'); % [loc,scale,nu] output
        for b=1:B
            param_bs_sample_matlab = mle_param_matlab(1)+mle_param_matlab(2)*trnd(mle_param_matlab(3),T_samples(t),1); 
            VaR_param_matlab = quantile(param_bs_sample_matlab, alpha); 
            temp = param_bs_sample_matlab(param_bs_sample_matlab<=VaR_param_matlab);
            ESvec_param_matlab(b) = mean(temp);
        end
        
        % length of CI (parametric bs)
        ci_param_matlab = quantile(ESvec_param_matlab, [alpha/2 1-alpha/2]);
        low_param_matlab = ci_param_matlab(1); high_param_matlab = ci_param_matlab(2);
        ci_length_param_matlab(r,t) = high_param_matlab-low_param_matlab;
        
        % whether or not the interval contains the TRUE ES (parametric bs)
        ci_trueES_param_matlab(r,t) = (trueES>low_param_matlab)&(trueES<high_param_matlab);
        %disp(ci_trueES_nonparam(r,t)
        
        %%% non-parametric bootstrap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % book code -> Listing 2.9, p.56
        for b=1:B
            ind = unidrnd(T_samples(t),[T_samples(t),1]);
            nonparam_bs_sample=data(ind);
            %phatvec_nonpara_bs(b) = 1/mean(nonparam_bs); % necessary?
            VaR_nonpara = quantile(nonparam_bs_sample, alpha); 
            temp_nonparam = nonparam_bs_sample(nonparam_bs_sample<=VaR_nonpara);
            ESvec_nonparam(b) = mean(temp_nonparam);
        end

        %%% length of CI (non-parametric bs)
        ci_nonparam = quantile(ESvec_nonparam, [alpha/2 1-alpha/2]);
        low_nonparam = ci_nonparam(1); high_nonparam = ci_nonparam(2);
        ci_length_nonparam(r,t) = high_nonparam-low_nonparam;
        
        %%% whether or not the interval contains the TRUE ES (non-parametric bs)
        ci_trueES_nonparam(r,t) = (trueES>low_nonparam)&(trueES<high_nonparam);
        %disp(ci_trueES_nonparam(r,t));
        
    end
    
    
    % Printing some stuff
    disp('================================')
    X = ['True ES:', num2str(trueES)]; 
    disp(X);

    text = ['For T = ', int2str(T_samples(t)),' Rep = ',int2str(rep),' B = ',int2str(B)];
    disp(text)
    disp('--------------------------------')
    disp('Parametric (Book codes):')
    mean_ci_length_param_book = mean(ci_length_param_book(:,t));
    A1 = ['Mean CI length:', num2str(mean_ci_length_param_book)];
    disp(A1)
    sum_bool1 = sum(ci_trueES_param_book(:,t));
    coverage_ratio1 = sum_bool1/length(ci_trueES_param_book(:,t));
    B1 = ['Coverage ratio:', num2str(coverage_ratio1)]; 
    disp(B1)

    disp('--------------------------------')
    disp('Parametric (MATLAB built-in):')
    mean_ci_length_param_matlab = mean(ci_length_param_matlab(:,t));
    A2 = ['Mean CI length:', num2str(mean_ci_length_param_matlab)];
    disp(A2)
    sum_bool2 = sum(ci_trueES_param_matlab(:,t));
    coverage_ratio2 = sum_bool2/length(ci_trueES_param_matlab(:,t));
    B2 = ['Coverage ratio:', num2str(coverage_ratio2)]; 
    disp(B2)

    disp('--------------------------------')
    disp('Non-parametric:')
    mean_ci_length_nonparam = mean(ci_length_nonparam(:,t));
    A3 = ['Mean CI length:', num2str(mean_ci_length_nonparam)]; 
    disp(A3)
    sum_bool3 = sum(ci_trueES_nonparam(:,t));
    coverage_ratio3 = sum_bool3/length(ci_trueES_nonparam(:,t));
    B3 = ['Coverage ratio:', num2str(coverage_ratio3)]; 
    disp(B3)


% additional check on the parameters 
    disp('--------------------------------')
    Z = ['BS parameters book MLE:', num2str(mle_param_book)]; 
    disp(Z)
    Y = ['BS parameters matlab MLE:', num2str(mle_param_matlab)]; 
    disp(Y)

end

%% Question 3 - Lo

% Repeat Q2 but using SYMMETRIC stable Paretian distribution. 
% Possible to use MATLAB's built in function for it. 
tic
% Parameters
a1 = 1.6; %a1 = 1.8; 
b = 0; 
c = 1;                      % scale
d = 0;                      % location
n = 1e7; 
xvec = -20:0.01:20; 
method = 1; 
T = 300;
loc = 0; scale = 1;         
rep = 300;                  % 1000 Repititions
alpha = 0.1; 
df = 3;                     % df = 6, df = 3;
T_samples = [100];            % [100 500 2000] Sample size
B = 300;                    % 1000 Bootstrap samples
seed = 2;                   % set seed value so we can replicate results
initvec = [df loc scale];   % for parametric MLE

%%% PART 1: simulation of a symmetric stable (meaning beta = 0) Paretian 
% distribution using alpha1 = 1.6 and alpha2 = 1.8. 


% Creating the kernel density estimate by generating 1e6 simulated IID
% variates from the stable distribution and simulating them. 
eststab = stabgen(n,a1,b,c,d,seed); 
[f,cd] = ksdensity(eststab,xvec); 
figure, hold on, plot(cd, f , 'r-', 'linewidth', 2)

t_dist = trnd(df,T,1);
[f2,cd2] = ksdensity(t_dist,xvec); 
plot(cd2, f2 , 'b--', 'linewidth', 2)

% adjusting the plot
xlim([-20 20])
legend('Simulated PDF of the Symmetric Stable Paretian distribution', 'Location', 'NorthWest')
title('Simulated PDF of the Symmetric Stable Paretian distribution')
hold off

%%% PART 2: compute the true ES of the symmetric stable distribution - 
% both via the Stoyanov et al. result and by simulation. 

% Theoretical ES based on a tail probability xi=0.01 using the Stoyanov et
% al. result
[Theo_ES, Theo_VaR] = asymstableES(alpha, a1, b, d, c, method); 
X = ['ES via Stoyanov et al: ', num2str(Theo_ES)];
disp(X);

% Simulation of ES
P = stabgen(n, a1, b, c, d, seed); 
VaR = quantile(P, alpha);
Plo=P(P<=VaR);
ES_simulated = mean(Plo);
X = ['Simulated ES: ' , num2str(ES_simulated)]; 
disp(X);

%%% PART 3: your codes for the parametric and nonparametric bootstrap 
% stay the same, but you USE AS THE TRUE DATA simulations from the 
% symmetric Stable Paretian distribution.


% Fixing the 3 true parameters: tailprob alpha, dof, left tail quantile of location 0, scale 1 distribution
%alpha = 0.1; df = 4; c01 = tinv(alpha , df); % alpha = 0.01


% Initializing vectors
ESvec_nonparam = zeros(B, 1); 
ESvec_param_book = zeros(B, 1);
ESvec_param_matlab = zeros(B, 1); 
ci_length_nonparam = zeros(rep, length(T_samples)); % length of ci interval for each rep
ci_length_param_book = zeros(rep, length(T_samples)); 
ci_length_param_matlab = zeros(rep, length(T_samples)); 
ci_trueES_nonparam = zeros(rep, length(T_samples)); % whether interval contains trueES for each rep
ci_trueES_param_book = zeros(rep, length(T_samples)); 
ci_trueES_param_matlab = zeros(rep, length(T_samples)); 


% True ES for student t
trueES = Theo_ES;
%trueES = ES_simulated;


% Simulate "rep" repetitions of an IID T-length sequence of Student t
% For each rep claculate bootstrap CI 90% based on B bootstrap replications
for t = 1:length(T_samples)
    for r = 1:rep
        %%% generating data points from student t distribution
        data = stabgen(T_samples(t),a1,b,c,d,seed);
        

        %%% parametric bootstrap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BOOK        
        % book code -> Listing 4.5, p.139
        mle_param_book = tlikmax0(data, initvec);

        % Now that I have the parameters' estimation through MLE based on the bootstrap data set
        % I want to compute the theoretical ES based on those parameter values. 
        for b=1:B
            param_bs_sample_book = mle_param_book(2)+mle_param_book(3)*trnd(mle_param_book(1),T_samples(t),1);
            VaR_param_book = quantile(param_bs_sample_book, alpha); 
            temp_book = param_bs_sample_book(param_bs_sample_book<=VaR_param_book); 
            ESvec_param_book(b) = mean(temp_book);
        end 

        % length of CI (parametric bs)
        ci_param_book = quantile(ESvec_param_book,[alpha/2 1-alpha/2]); 
        %ci_param_book = quantile(phatvec, [alpha/2 1-alpha/2]);
        low_param_book = ci_param_book(1); high_param_book = ci_param_book(2); 
        ci_length_param_book(r,t) = high_param_book-low_param_book; 

        % whether or not the interval contains the TRUE ES (parametric bs)
        ci_trueES_param_book(r,t) = (trueES>low_param_book)&(trueES<high_param_book);
        
        % MATLAB BUILT IN:
        mle_param_matlab = mle(data, 'Distribution', 'tLocationScale'); % [loc,scale,nu] output
        for b=1:B
            param_bs_sample_matlab = mle_param_matlab(1)+mle_param_matlab(2)*trnd(mle_param_matlab(3),T_samples(t),1); 
            VaR_param_matlab = quantile(param_bs_sample_matlab, alpha); 
            temp = param_bs_sample_matlab(param_bs_sample_matlab<=VaR_param_matlab);
            ESvec_param_matlab(b) = mean(temp);
        end
        
        % length of CI (parametric bs)
        ci_param_matlab = quantile(ESvec_param_matlab, [alpha/2 1-alpha/2]);
        low_param_matlab = ci_param_matlab(1); high_param_matlab = ci_param_matlab(2);
        ci_length_param_matlab(r,t) = high_param_matlab-low_param_matlab;
        
        % whether or not the interval contains the TRUE ES (parametric bs)
        ci_trueES_param_matlab(r,t) = (trueES>low_param_matlab)&(trueES<high_param_matlab);
        %disp(ci_trueES_nonparam(r,t)
        
        %%% non-parametric bootstrap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % book code -> Listing 2.9, p.56
        for b=1:B
            ind = unidrnd(T_samples(t),[T_samples(t),1]);
            nonparam_bs_sample=data(ind);
            %phatvec_nonpara_bs(b) = 1/mean(nonparam_bs); % necessary?
            VaR_nonpara = quantile(nonparam_bs_sample, alpha); 
            temp_nonparam = nonparam_bs_sample(nonparam_bs_sample<=VaR_nonpara);
            ESvec_nonparam(b) = mean(temp_nonparam);
        end

        %%% length of CI (non-parametric bs)
        ci_nonparam = quantile(ESvec_nonparam, [alpha/2 1-alpha/2]);
        low_nonparam = ci_nonparam(1); high_nonparam = ci_nonparam(2);
        ci_length_nonparam(r,t) = high_nonparam-low_nonparam;
        
        %%% whether or not the interval contains the TRUE ES (non-parametric bs)
        ci_trueES_nonparam(r,t) = (trueES>low_nonparam)&(trueES<high_nonparam);
        %disp(ci_trueES_nonparam(r,t));
        
    end
    
    
    % Printing some stuff
    disp('================================')
    X = ['True ES:', num2str(trueES)]; 
    disp(X);

    text = ['For T = ', int2str(T_samples(t)),' Rep = ',int2str(rep),' B = ',int2str(B)];
    disp(text)
    disp('--------------------------------')
    disp('Parametric (Book codes):')
    mean_ci_length_param_book = mean(ci_length_param_book(:,t));
    A1 = ['Mean CI length:', num2str(mean_ci_length_param_book)];
    disp(A1)
    sum_bool1 = sum(ci_trueES_param_book(:,t));
    coverage_ratio1 = sum_bool1/length(ci_trueES_param_book(:,t));
    B1 = ['Coverage ratio:', num2str(coverage_ratio1)]; 
    disp(B1)

    disp('--------------------------------')
    disp('Parametric (MATLAB built-in):')
    mean_ci_length_param_matlab = mean(ci_length_param_matlab(:,t));
    A2 = ['Mean CI length:', num2str(mean_ci_length_param_matlab)];
    disp(A2)
    sum_bool2 = sum(ci_trueES_param_matlab(:,t));
    coverage_ratio2 = sum_bool2/length(ci_trueES_param_matlab(:,t));
    B2 = ['Coverage ratio:', num2str(coverage_ratio2)]; 
    disp(B2)

    disp('--------------------------------')
    disp('Non-parametric:')
    mean_ci_length_nonparam = mean(ci_length_nonparam(:,t));
    A3 = ['Mean CI length:', num2str(mean_ci_length_nonparam)]; 
    disp(A3)
    sum_bool3 = sum(ci_trueES_nonparam(:,t));
    coverage_ratio3 = sum_bool3/length(ci_trueES_nonparam(:,t));
    B3 = ['Coverage ratio:', num2str(coverage_ratio3)]; 
    disp(B3)


% additional check on the parameters 
    disp('--------------------------------')
    Z = ['BS parameters book MLE:', num2str(mle_param_book)]; 
    disp(Z)
    Y = ['BS parameters matlab MLE:', num2str(mle_param_matlab)]; 
    disp(Y)

end

toc
%% Question 4:

%%%% Idea: basically repeat Question 1, but using the NCT for both the true DGP, *and also* use the NCT for the parametric bootstrap
% -> use same parameter values as in Question 2: 
    % dof 3 and 6, 3 sample size T, set of noncentrality parameter mu values (I had suggested -3, -2, -1, 0). 

% compute pdf of the location-zero, scale-one NCT -> Fundamental Statistics, section 9.3.2,

% compute the MLE of the location-scale NCT -> Fundamental Statistics, chapter  (two ways)

    % via Paolella's "d.d.a" NCT approximation 
    
    % via MATLAB's built in pdf of NCT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
mu = -3; % [0 -1 -2 -3]  
df = 3; % [3 6]         
x = (-20:0.1:20)'; 

alpha = 0.1; 
n = 10^6;

loc = 0; scale = 1;         % loc = 1; scale = 2;
rep = 1000;                  % 1000 Repititions
T_samples = [500];           % [100 500 2000] Sample size
B = 1000;                    % 1000 Bootstrap samples
seed = 2;         % set seed value so we can replicate results
initvec = [df loc scale mu];   % for parametric MLE
method = 'matlab';             % ddanct or matlab


% compute pdf of the location-zero, scale-one NCT -> Fundamental Statistics, section 9.3.2, 
%nct_book = stdnctpdfln_j(x, df, mu); 
nct_book = exp(stdnctpdfln_j(x, df, mu)); 

% MATLAB built in function --> works, but we have to compare both 
% (actually we don't)
%nct_matlab = log(nctpdf(x, df, mu)); 
nct_matlab = nctpdf(x, df, mu); 

figure, plot(x, nct_book, 'b','LineWidth', 2)
hold on, plot(x, nct_matlab, 'r--','LineWidth', 2)
xlim([-20 20])
legend('Book code','MATLAB built-in')
title('Simulated PDF of the NCT')
hold off




% ES via simulation
%nct_rnd_sample = loc+scale*nctrnd(df,mu,100,1);
norm=normrnd(mu,1,n,1); chi2=chi2rnd(df,n,1);         % central chi2
nct_rnd_sample = loc+scale*(norm./sqrt(chi2/df));
VaR = quantile(nct_rnd_sample, alpha);
temp_2 = nct_rnd_sample(nct_rnd_sample <= VaR);
ES_simulated = mean(temp_2); 
X = ['Simulated ES: ', num2str(ES_simulated)]; 
disp(X);

% ES using integral definition of NCT (Ila)
c01_nct = nctinv(alpha , df, mu); 
ES_01_nct = @(x) x.*nctpdf(x, df, mu); 
ES_nct_int= integral(ES_01_nct, -Inf, c01_nct)/alpha; 
disp(['ES via Numeric Integration:', num2str(ES_nct_int)]);

% True ES for student t
%trueES = ES_simulated; 
trueES = ES_nct_int;


% Initializing vectors
ESvec_nonparam = zeros(B, 1); 
ESvec_param_book = zeros(B, 1); 
ESvec_param_matlab = zeros(B, 1); 
ci_length_nonparam = zeros(rep, length(T_samples)); % length of ci interval for each rep
ci_length_param_book = zeros(rep, length(T_samples)); 
ci_length_param_matlab = zeros(rep, length(T_samples)); 
ci_trueES_nonparam = zeros(rep, length(T_samples)); % whether interval contains trueES for each rep
ci_trueES_param_book = zeros(rep, length(T_samples)); 
ci_trueES_param_matlab = zeros(rep, length(T_samples)); 


% Simulate "rep" repetitions of an IID T-length sequence of Student t
% For each rep claculate bootstrap CI 90% based on B bootstrap replications
for t = 1:length(T_samples)
    for r = 1:rep
        %%% generating data points from student t distribution    
        norm=normrnd(mu,1, T_samples(t), 1); chi2=chi2rnd(df, T_samples(t), 1);  % noncentral chi2
        data = loc+scale*(norm./sqrt(chi2/df));
        
        
        %%% parametric bootstrap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BOOK
        % book code -> use Listing 4.5/4.6 p.141 
        % method = ddnct or matlab ("d.d.a" NCT & Matlab built-in)
        mle_param_book = tlikmax0_modified(data, initvec, method);

        % Now that I have the parameters' estimation through MLE based on the bootstrap data set
        % I want to compute the theoretical ES based on those parameter values. 

        for b=1:B
            rnd_norm = normrnd(mle_param_book(4),1, T_samples(t), 1); rnd_chi2=chi2rnd(mle_param_book(1), T_samples(t), 1);
            rnd_nct = rnd_norm./sqrt(rnd_chi2/mle_param_book(1));
            param_bs_sample_book = mle_param_book(2)+mle_param_book(3)*rnd_nct;
            VaR_param_book = quantile(param_bs_sample_book, alpha); 
            temp_book = param_bs_sample_book(param_bs_sample_book<=VaR_param_book); 
            ESvec_param_book(b) = mean(temp_book);
        end 

        % length of CI (parametric bs)
        ci_param_book = quantile(ESvec_param_book,[alpha/2 1-alpha/2]); 
        %ci_param_book = quantile(phatvec, [alpha/2 1-alpha/2]);
        low_param_book = ci_param_book(1); high_param_book = ci_param_book(2); 
        ci_length_param_book(r,t) = high_param_book-low_param_book; 

        % whether or not the interval contains the TRUE ES (parametric bs)
        ci_trueES_param_book(r,t) = (trueES>low_param_book)&(trueES<high_param_book); 

        
        
%         % MATLAB BUILT IN:
%         mle_param_matlab = mle(data, 'Distribution', 'tLocationScale'); % [loc,scale,nu] output
%         for b=1:B
%             param_bs_sample_matlab = mle_param_matlab(1)+mle_param_matlab(2)*trnd(mle_param_matlab(3),T_samples(t),1); 
%             VaR_param_matlab = quantile(param_bs_sample_matlab, alpha); 
%             temp = param_bs_sample_matlab(param_bs_sample_matlab<=VaR_param_matlab);
%             ESvec_param_matlab(b) = mean(temp);
%         end
%         
%         % length of CI (parametric bs)
%         ci_param_matlab = quantile(ESvec_param_matlab, [alpha/2 1-alpha/2]);
%         low_param_matlab = ci_param_matlab(1); high_param_matlab = ci_param_matlab(2);
%         ci_length_param_matlab(r,t) = high_param_matlab-low_param_matlab;
%         
%         % whether or not the interval contains the TRUE ES (parametric bs)
%         ci_trueES_param_matlab(r,t) = (trueES>low_param_matlab)&(trueES<high_param_matlab);
%         
        
        
        
        %%% non-parametric bootstrap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % book code -> Listing 2.9, p.56
        for b=1:B
            ind = unidrnd(T_samples(t),[T_samples(t),1]);
            nonparam_bs_sample=data(ind);
            %phatvec_nonpara_bs(b) = 1/mean(nonparam_bs); % necessary?
            VaR_nonpara = quantile(nonparam_bs_sample, alpha); 
            temp_nonparam = nonparam_bs_sample(nonparam_bs_sample<=VaR_nonpara);
            ESvec_nonparam(b) = mean(temp_nonparam);
        end

        %%% length of CI (non-parametric bs)
        ci_nonparam = quantile(ESvec_nonparam, [alpha/2 1-alpha/2]);
        low_nonparam = ci_nonparam(1); high_nonparam = ci_nonparam(2);
        ci_length_nonparam(r,t) = high_nonparam-low_nonparam;
        
        %%% whether or not the interval contains the TRUE ES (non-parametric bs)
        ci_trueES_nonparam(r,t) = (trueES>low_nonparam)&(trueES<high_nonparam);
        %disp(ci_trueES_nonparam(r,t));
        
        
    end
    
    
    % Printing some stuff
    disp('================================')
    X = ['True ES:', num2str(trueES)]; 
    disp(X);

    text = ['For T = ', int2str(T_samples(t)),' Rep = ',int2str(rep),' B = ',int2str(B)];
    disp(text)
    disp('--------------------------------')
    disp('Parametric (Book codes):')
    mean_ci_length_param_book = mean(ci_length_param_book(:,t));
    A1 = ['Mean CI length:', num2str(mean_ci_length_param_book)];
    disp(A1)
    sum_bool1 = sum(ci_trueES_param_book(:,t));
    coverage_ratio1 = sum_bool1/length(ci_trueES_param_book(:,t));
    B1 = ['Coverage ratio:', num2str(coverage_ratio1)]; 
    disp(B1)

%     disp('--------------------------------')
%     disp('Parametric (MATLAB built-in):')
%     mean_ci_length_param_matlab = mean(ci_length_param_matlab(:,t));
%     A2 = ['Mean CI length:', num2str(mean_ci_length_param_matlab)];
%     disp(A2)
%     sum_bool2 = sum(ci_trueES_param_matlab(:,t));
%     coverage_ratio2 = sum_bool2/length(ci_trueES_param_matlab(:,t));
%     B2 = ['Coverage ratio:', num2str(coverage_ratio2)]; 
%     disp(B2)

    disp('--------------------------------')
    disp('Non-parametric:')
    mean_ci_length_nonparam = mean(ci_length_nonparam(:,t));
    A3 = ['Mean CI length:', num2str(mean_ci_length_nonparam)]; 
    disp(A3)
    sum_bool3 = sum(ci_trueES_nonparam(:,t));
    coverage_ratio3 = sum_bool3/length(ci_trueES_nonparam(:,t));
    B3 = ['Coverage ratio:', num2str(coverage_ratio3)]; 
    disp(B3)
    
    % additional check on the parameters 
    disp('--------------------------------')
    Z = ['BS parameters book MLE:', num2str(mle_param_book)]; 
    disp(Z)
%     Y = ['BS parameters matlab MLE:', num2str(mle_param_matlab)]; 
%     disp(Y)



end
