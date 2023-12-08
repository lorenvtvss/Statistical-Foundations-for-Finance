function [final_nu, nu_vec, mu, sigma] = function_MMFAlgorithm(x, initial_df, weights, reps)
% input parameters
% % x                   matrix of random samples of a multivarite t dist
% % initial_df          starting value for the degrees of freedom
% % weights             weights
% % reps                number of repetitions
%
% output
% % nu_vec              estimate of the degrees of freedom of a Student t distribution
% % mu                  mean vector of the latest iteration
% % sigma               variance-covariance matrix of the latest iteration


% check input
if initial_df <= 0
    error ("df_initial must be strictly larger 0")
end

if size(x, 1) > size(x, 2)
    error("number of samples must be less than dim + 1 (where dim: sample size)")
end

if sum(weights <= 0) > 0
    error("all weights must be strictly positive")
end

if sum(weights) - 1 > 1e-10
    error("the weights must sum to one")
end

if sum(weights) - 1 < -1e-10
    error("the weights must sum to one")
end



% initialize variables
d = size(x , 1);            % dimension=3
nu = initial_df;            % nu
nu_vec = zeros(reps, 1);
mu = sum(x, 2)/size(x, 2);  % mu

% initialize sigma matrix
sigma0 = 0;
for i = 1:length(x)
    sigma0 = sigma0 + (x(:, i) - mu)*(x(:, i) - mu)';
end
sigma = sigma0/length(x);


for r = 1:reps
    %%%%% E-step: compute weights
    % initialize vectors for delta and gamma
    delta = zeros(size(x, 2), 1);
    gamma = zeros(size(x, 2), 1);
    % fill vectors with values for the current rep loop
    for i = 1:size(x, 2)
       delta(i) = (x(:, i) - mu)' / sigma * (x(:, i) - mu);
       gamma(i) = (nu + d)/(nu + delta(i,1));
    end

    %%%%% M-step: update the parameters
    % initialize (set to zero) variable to save denominator for updating mu and sigma
    denom = 0;
    for i = 1:size(x, 2)
        denom = denom + weights(i)*gamma(i);
    end

    %%% mu
    % initialize (set to zero) variable to save the numerator
    mu_nom = 0;
    % calculate the nominator for updating mu
    for i = 1:size(x, 2)
        mu_nom = mu_nom + weights(i) * gamma(i) * x(:, i);
    end
    % update mu value
    mu = mu_nom / denom;

    %%% sigma
    % initialize (set to zero) variable to save the numerator
    sigma_num = 0;
    % calculate the nominator for updating sigma
    for i = 1:length(x)
        sigma_num = sigma_num + ( weights(i) * gamma(i) * (x(:,i) - mu) * (x(:,i) - mu)' );
    end
    % update sigma value
    sigma = sigma_num / denom;

    %%% nu
    % initialize (set to zero) variable to save the sum part of the updating step for nu
    nu_sum = 0;
    for i = 1:length(x)
        nu_sum = nu_sum + weights(i) * ( (nu + d) / (nu + delta(i)) - log( (nu + d)/(nu + delta(i)) ) - 1 );
    end
    nu = fzero(@(x)  phi_func(x/2) - phi_func((x+d) / 2) + nu_sum , [1e-100, 1e100]);
    nu_vec(r) = nu;

end 
final_nu = nu;

end 

function [phi] = phi_func(x)
    % see p. 81
    phi = psi(x) - log(x);
end