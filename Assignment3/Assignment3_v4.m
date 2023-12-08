%% Assignment 1 - complete

%%% Create covariance matrix with positive eigenvalues and diagonal values
%%% of 1
% Input parameters        
dim = 3;            % matrix dimensionality
diag = 1;           % diagonal values
M = zeros(3, 3);

while (sum(eig(M)<=0) ~= 0)
    t = triu(bsxfun(@min,diag,-diag.').*rand(dim),1); % The upper trianglar random values
    M = diag(diag)+t+t.'; % Put them together in a symmetric matrix
end
disp('Sigma matrix: '); disp(M);
disp('Eigenvalues of sigma matrix: '); disp(eig(M));

\text
%%% Simulate 3d MVT and estimate parameters

% Parameters
T_samples = [200,2000];    % [200, 2000]
rep = 500;                 % # repititions: 500
loc = 0;
scale = 1;
df_true = 4;               % true df value
reps = 500;                % # repititions for MMF: 500
initial_df = 0.1;

% Parameters for MMF
step_algorithm = 'MMF';
anz_steps = 500;
stop = 1;
abs_criteria = 1;
regularize = 1;
save_obj = 1;


% Initializing vectors
mus_mmf = zeros(rep, length(T_samples)*dim); 
nus_mmf = zeros(rep, length(T_samples)); 
sigmas_mmf = zeros(rep, length(T_samples)*6); 
mus_mmf_adj = zeros(rep, length(T_samples)*dim); 
nus_mmf_adj = zeros(rep, length(T_samples)); 
sigmas_mmf_adj = zeros(rep, length(T_samples)*6);

mus_mle = zeros(rep, length(T_samples)*dim); 
nus_mle = zeros(rep, length(T_samples)); 
sigmas_mle = zeros(rep, length(T_samples)*6); 
mus_mle_adj = zeros(rep, length(T_samples)*dim); 
nus_mle_adj = zeros(rep, length(T_samples)); 
sigmas_mle_adj = zeros(rep, length(T_samples)*6);

mus_ecme = zeros(rep, length(T_samples)*dim); 
nus_ecme = zeros(rep, length(T_samples)); 
sigmas_ecme = zeros(rep, length(T_samples)*6); 
mus_ecme_adj = zeros(rep, length(T_samples)*dim); 
nus_ecme_adj = zeros(rep, length(T_samples)); 
sigmas_ecme_adj = zeros(rep, length(T_samples)*6);

mus_approx = zeros(rep, length(T_samples)*dim); 
nus_approx = zeros(rep, length(T_samples)); 
sigmas_approx = zeros(rep, length(T_samples)*6); 
mus_approx_adj = zeros(rep, length(T_samples)*dim); 
nus_approx_adj = zeros(rep, length(T_samples)); 
sigmas_approx_adj = zeros(rep, length(T_samples)*6);

time_mmf = zeros(1,length(T_samples))
time_mle = zeros(1,length(T_samples))
time_ecme = zeros(1,length(T_samples))
time_approx = zeros(1,length(T_samples))

ts = 0; tm = 0;

% True values
mu_true = [0 0 0]; 
df_true = 4;
sigma_true = M;

for t = 1:length(T_samples)
    
    % Simulate multi variate t distribution (MVT)
    data=loc+scale*mvtrnd(M,df_true,T_samples(t)); 
    
    for r = 1:rep
        
        % Sampling from MVT
        ind = unidrnd(T_samples(t),[T_samples(t),1]);
        mvt_sample=data(ind,1:3)';

        %%% Estimating parameters through MMF
        % and creating adjusted estimates by true values
        weights = 1/T_samples(t) * ones(T_samples(t), 1);
        tic
        [nu, nu_vec, mu, sigma] = function_MMFAlgorithm(mvt_sample, initial_df, weights, reps);
        time_mmf(t) = toc;
        
        mus_mmf(r,1+tm) = mu(1);   mus_mmf_adj(r,1+tm) = mu(1) - mu_true(1);
        mus_mmf(r,2+tm) = mu(2);   mus_mmf_adj(r,2+tm) = mu(2) - mu_true(2);
        mus_mmf(r,3+tm) = mu(3);   mus_mmf_adj(r,3+tm) = mu(3) - mu_true(3);
        
        nus_mmf(r,t) = nu;   nus_mmf_adj(r,t) = nu - df_true;
        
        sigmas_mmf(r,1+ts) = sigma(1,1);  sigmas_mmf_adj(r,1+ts) = sigma(1,1) - M(1,1);
        sigmas_mmf(r,2+ts) = sigma(2,1);  sigmas_mmf_adj(r,2+ts) = sigma(2,1) - M(2,1);
        sigmas_mmf(r,3+ts) = sigma(3,1);  sigmas_mmf_adj(r,3+ts) = sigma(3,1) - M(3,1);
        sigmas_mmf(r,4+ts) = sigma(2,2);  sigmas_mmf_adj(r,4+ts) = sigma(2,2) - M(2,2);
        sigmas_mmf(r,5+ts) = sigma(3,2);  sigmas_mmf_adj(r,5+ts) = sigma(3,2) - M(3,2);
        sigmas_mmf(r,6+ts) = sigma(3,3);  sigmas_mmf_adj(r,6+ts) = sigma(3,3) - M(3,3);
        
        %%% Estimating parameters through MLE
        % and creating adjusted estimates by true values
        tic
        [param,stderr,iters,loglik,Varcov] = MVTestimation_3d(mvt_sample', weights);
        time_mle(t) = toc;
        
        % param output order: 
        % k, mu1, mu2, mu3, Sigma_11, Sigma_12, Sigma_13,Sigma_22, Sigma_23, Sigma_33)
        mus_mle(r,1+tm) = param(2);   mus_mle_adj(r,1+tm) = param(2) - mu_true(1);
        mus_mle(r,2+tm) = param(3);   mus_mle_adj(r,2+tm) = param(3) - mu_true(2);
        mus_mle(r,3+tm) = param(4);   mus_mle_adj(r,3+tm) = param(4) - mu_true(3);
        
        nus_mle(r,t) = param(1);   nus_mle_adj(r,t) = param(1) - df_true;
        
        sigmas_mle(r,1+ts) = param(5);  sigmas_mle_adj(r,1+ts) = param(5) - M(1,1);
        sigmas_mle(r,2+ts) = param(6);  sigmas_mle_adj(r,2+ts) = param(6) - M(2,1);
        sigmas_mle(r,3+ts) = param(7);  sigmas_mle_adj(r,3+ts) = param(7) - M(3,1);
        sigmas_mle(r,4+ts) = param(8);  sigmas_mle_adj(r,4+ts) = param(8) - M(2,2);
        sigmas_mle(r,5+ts) = param(9);  sigmas_mle_adj(r,5+ts) = param(9) - M(3,2);
        sigmas_mle(r,6+ts) = param(10);  sigmas_mle_adj(r,6+ts) = param(10) - M(3,3);

        %%% Estimating parameters through ECME
        % and creating adjusted estimates by true values
        tic
        [mu_ECME, S_ECME, nu_ECME] = fitt(mvt_sample');
        time_ecme(t) = toc;
        
        mus_ecme(r,1+tm) = mu_ECME(1);   mus_ecme_adj(r,1+tm) = mu_ECME(1) - mu_true(1);
        mus_ecme(r,2+tm) = mu_ECME(2);   mus_ecme_adj(r,2+tm) = mu_ECME(2) - mu_true(2);
        mus_ecme(r,3+tm) = mu_ECME(3);   mus_ecme_adj(r,3+tm) = mu_ECME(3) - mu_true(3);
        
        nus_ecme(r,t) = nu_ECME;   nus_ecme_adj(r,t) = nu_ECME - df_true;
        
        sigmas_ecme(r,1+ts) = S_ECME(1,1);  sigmas_ecme_adj(r,1+ts) = S_ECME(1,1) - M(1,1);
        sigmas_ecme(r,2+ts) = S_ECME(2,1);  sigmas_ecme_adj(r,2+ts) = S_ECME(2,1) - M(2,1);
        sigmas_ecme(r,3+ts) = S_ECME(3,1);  sigmas_ecme_adj(r,3+ts) = S_ECME(3,1) - M(3,1);
        sigmas_ecme(r,4+ts) = S_ECME(2,2);  sigmas_ecme_adj(r,4+ts) = S_ECME(2,2) - M(2,2);
        sigmas_ecme(r,5+ts) = S_ECME(3,2);  sigmas_ecme_adj(r,5+ts) = S_ECME(3,2) - M(3,2);
        sigmas_ecme(r,6+ts) = S_ECME(3,3);  sigmas_ecme_adj(r,6+ts) = S_ECME(3,3) - M(3,3);

        %%% Estimating parameters through the approximation method 
        % and creating adjusted estimates by true values
        tic
        [mu_approx, S_approx, nu_approx] = fitt_approx(mvt_sample'); 
        time_approx(t) = toc;

        mus_approx(r,1+tm) = mu_approx(1);   mus_approx_adj(r,1+tm) = mu_approx(1) - mu_true(1);
        mus_approx(r,2+tm) = mu_approx(2);   mus_approx_adj(r,2+tm) = mu_approx(2) - mu_true(2);
        mus_approx(r,3+tm) = mu_approx(3);   mus_approx_adj(r,3+tm) = mu_approx(3) - mu_true(3);
        
        nus_approx(r,t) = nu_approx;   nus_approx_adj(r,t) = nu_approx - df_true;
        
        sigmas_approx(r,1+ts) = S_approx(1,1);  sigmas_approx_adj(r,1+ts) = S_approx(1,1) - M(1,1);
        sigmas_approx(r,2+ts) = S_approx(2,1);  sigmas_approx_adj(r,2+ts) = S_approx(2,1) - M(2,1);
        sigmas_approx(r,3+ts) = S_approx(3,1);  sigmas_approx_adj(r,3+ts) = S_approx(3,1) - M(3,1);
        sigmas_approx(r,4+ts) = S_approx(2,2);  sigmas_approx_adj(r,4+ts) = S_approx(2,2) - M(2,2);
        sigmas_approx(r,5+ts) = S_approx(3,2);  sigmas_approx_adj(r,5+ts) = S_approx(3,2) - M(3,2);
        sigmas_approx(r,6+ts) = S_approx(3,3);  sigmas_approx_adj(r,6+ts) = S_approx(3,3) - M(3,3);
    end
    ts = ts+6; tm = tm+dim;
end
%% Printing values
disp('********MMF********')
disp(['Sample size: ', num2str(T_samples(1))])
fprintf('Estimation time: %d min %f sec\n', floor(time_mmf(1)/60), rem(time_mmf(1),60));
disp(['Mean nu: ', num2str(mean(nus_mmf(:,1)))])
disp(['Mean mu: ', num2str(mean(mus_mmf(:,1:3)))])
disp(['Sample size: ', num2str(T_samples(2))])
fprintf('Estimation time: %d min %f sec\n', floor(time_mmf(2)/60), rem(time_mmf(2),60));
disp(['Mean nu: ', num2str(mean(nus_mmf(:,2)))])
disp(['Mean mu: ', num2str(mean(mus_mmf(:,4:6)))])

disp('********MLE********')
disp(['Sample size: ', num2str(T_samples(1))])
fprintf('Estimation time: %d min %f sec\n', floor(time_mle(1)/60), rem(time_mle(1),60));
disp(['Mean nu: ', num2str(mean(nus_mle(:,1)))])
disp(['Mean mu: ', num2str(mean(mus_mle(:,1:3)))])
disp(['Sample size: ', num2str(T_samples(2))])
fprintf('Estimation time: %d min %f sec\n', floor(time_mle(2)/60), rem(time_mle(2),60));
disp(['Mean nu: ', num2str(mean(nus_mle(:,2)))])
disp(['Mean mu: ', num2str(mean(mus_mle(:,4:6)))])

disp('*****ECME******') 
disp(['Sample size: ', num2str(T_samples(1))])
fprintf('Estimation time: %d min %f sec\n', floor(time_ecme(1)/60), rem(time_ecme(1),60));
disp(['Mean nu: ', num2str(mean(nus_ecme(:,1)))])
disp(['Mean mu: ', num2str(mean(mus_ecme(:,1:3)))])
disp(['Sample size: ', num2str(T_samples(2))])
fprintf('Estimation time: %d min %f sec\n', floor(time_ecme(2)/60), rem(time_ecme(2),60));
disp(['Mean nu: ', num2str(mean(nus_ecme(:,2)))])
disp(['Mean mu: ', num2str(mean(mus_ecme(:,4:6)))])

disp('*****APPROX method******') 
disp(['Sample size: ', num2str(T_samples(1))])
fprintf('Estimation time: %d min %f sec\n', floor(time_approx(1)/60), rem(time_approx(1),60));
disp(['Mean nu: ', num2str(mean(nus_approx(:,1)))])
disp(['Mean mu: ', num2str(mean(mus_approx(:,1:3)))])
disp(['Sample size: ', num2str(T_samples(2))])
fprintf('Estimation time: %d min %f sec\n', floor(time_approx(2)/60), rem(time_approx(2),60));
disp(['Mean nu: ', num2str(mean(nus_approx(:,2)))])
disp(['Mean mu: ', num2str(mean(mus_approx(:,4:6)))])




%% Grouped plot - mus
hAxes.TickLabelInterpreter = 'latex';

figure, lab = {'\mu_1', '\mu_2', '\mu_3'}; 
boxplotGroup({mus_mmf_adj(:,1:3),mus_mle_adj(:,1:3),mus_ecme_adj(:,1:3),mus_approx_adj(:,1:3)}, 'PrimaryLabels',{'MMF','MLE','ECME','Aeschlimann'},'SecondaryLabels',lab,'GroupLabelType','Vertical','Whisker',1.5); set(gca,'fontsize',12)
title(['\mu values (sample size: ',num2str(T_samples(1)),')'])
name_mu_1 = ['Assignment3_ex1_mu_',T_samples(1),'.png'];
saveas(gcf,name_mu_1)

figure, lab = {'\mu_1', '\mu_2', '\mu_3'}; 
boxplotGroup({mus_mmf_adj(:,4:6),mus_mle_adj(:,4:6),mus_ecme_adj(:,4:6),mus_approx_adj(:,4:6)}, 'PrimaryLabels',{'MMF','MLE','ECME','Aeschlimann'},'SecondaryLabels',lab,'GroupLabelType','Vertical','Whisker',1.5); set(gca,'fontsize',12)
title(['\mu values (sample size: ',num2str(T_samples(2)),')'])
name_mu_2 = ['Assignment3_ex1_mu_',T_samples(2),'.png'];
saveas(gcf,name_mu_2)


%% Grouped plot - sigma
figure, lab = {'\sigma_{11}', '\sigma_{21}','\sigma_{31}', '\sigma_{22}','\sigma_{32}', '\sigma_{33}'}; 
boxplotGroup({sigmas_mmf_adj(:,1:6),sigmas_mle_adj(:,1:6),sigmas_ecme_adj(:,1:6),sigmas_approx_adj(:,1:6)},'PrimaryLabels',{'MMF','MLE','ECME','Aeschlimann'},'SecondaryLabels',lab,'GroupLabelType','Vertical','Whisker',1.5);
title(['\Sigma values (sample size: ',num2str(T_samples(1)),')'])
name_sigma_1 = ['Assignment3_ex1_sigma_',T_samples(1),'.png'];
saveas(gcf,name_sigma_1)


figure, lab = {'\sigma_{11}', '\sigma_{21}','\sigma_{31}', '\sigma_{22}','\sigma_{32}', '\sigma_{33}'}; 
boxplotGroup({sigmas_mmf_adj(:,7:12),sigmas_mle_adj(:,7:12),sigmas_ecme_adj(:,7:12),sigmas_approx_adj(:,7:12)},'PrimaryLabels',{'MMF','MLE','ECME','Aeschlimann'},'SecondaryLabels',lab,'GroupLabelType','Vertical','Whisker',1.5);
title(['\Sigma values (sample size: ',num2str(T_samples(2)),')'])
name_sigma_2 = ['Assignment3_ex1_sigma_',T_samples(2),'.png'];
saveas(gcf,name_sigma_2)


%% Grouped plot - nus
% not adjusted values
figure
boxplot([nus_mmf(:,1),nus_mle(:,1),nus_ecme(:,1),nus_approx(:,1)],'Labels',{'MMF','MLE','ECME','Aeschlimann'},'Whisker',1.5)
set(gca,'fontsize',12)
title(['\nu values (sample size: ',num2str(T_samples(1)),')'])
name_nu_1 = ['Assignment3_ex1_nu_',T_samples(1),'.png'];
saveas(gcf,name_nu_1)

figure
boxplot([nus_mmf(:,2),nus_mle(:,2),nus_ecme(:,2),nus_approx(:,2)],'Labels',{'MMF','MLE','ECME','Aeschlimann'},'Whisker',1.5)
set(gca,'fontsize',12)
title(['\nu values (sample size: ',num2str(T_samples(2)),')'])
name_nu_2 = ['Assignment3_ex1_nu_',T_samples(2),'.png'];
saveas(gcf,name_nu_2)





