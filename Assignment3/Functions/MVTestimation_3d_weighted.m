function [param,stderr,iters,loglik,Varcov] = MVTestimation_3d_weighted(x, rho)

% Program Listing 13.1 - Linear Models and Time series
% compute the weight vector
T = length(x); tvec = (1:T); omega=(T-tvec + 1).^(rho-1); weights = omega'/sum(omega);
disp(['Sum of weights: ', num2str(sum(weights))]);

% we recall the original MVT algorithm
[param,stderr,iters,loglik,Varcov] = MVTestimation_3d(x, weights)

end

