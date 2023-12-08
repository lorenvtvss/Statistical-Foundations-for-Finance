function [final_nu, nu_vec, mu, sigma] = MMFAlgorithm_weighted(rho, x, initial_df, reps)

% Program Listing 13.1 - Linear Models and Time series
% compute the weight vector
T = length(x); tvec = (1:T); omega=(T-tvec + 1).^(rho-1); weights = omega'/sum(omega);
disp(['Sum of weights: ', num2str(sum(weights))]);

% we recall the original MMF algorithm
[final_nu, nu_vec, mu, sigma] = function_MMFAlgorithm(x, initial_df, weights, reps);
end
