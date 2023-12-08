mu = 0; % [0 -1 -2 -3]  
df = 3; % [3 6]         
x = (-20:0.1:20)'; 

alpha = 0.1; 
n = 10^6;

loc = 0; scale = 1;         % loc = 1; scale = 2;
rep = 100;                  % 1000 Repititions
T = 100;           % [100 500 2000] Sample size
B = 100;                    % 1000 Bootstrap samples
seed = 2;         % set seed value so we can replicate results
%initvec = [df loc scale mu];   % for parametric MLE
initvec = [df loc scale ];
method = 'matlab';  
   

norm=normrnd(mu,1, T, 1); chi2=chi2rnd(df, T, 1);  % noncentral chi2
data = loc+scale*(norm./sqrt(chi2/df));

%mle_param_book = tlikmax0_modified(data, initvec, method);
mle_param_book = tlikmax0(data, initvec);