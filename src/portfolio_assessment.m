%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get asset return data
%
% INPUT
%   T1: Estimation window size
%   smonth: out-of-sample start month
%   emonth: out-of-sample end month
% OUTPUT
%   r: (T1+T2)xN asset returns (in excess of risk-free rate)
%   rf: (T1+T2)x1 risk-free returns
%   r2: T2xN out-of-sample asset returns (in excess of risk-free rate)
%   rf2: T2x1 out-of-sample risk-free returns
%   * T1, T2: estimation window size and out-of-sample period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [avg_utility ,avg_return, std_deviation ,sharpe_ratio] = portfolio_assessment(T1,gamma,with_transaction_cost)

% 0. define start/end date and other parameters 
smonth = 197601;
emonth = 201512;
transaction_cost_ratio = 0.003;
%T1 = 120;
%gamma = 1; % relative risk aversion coefficient

[r, rf, r2, rf2] = GetData(T1, smonth, emonth);
[T2,N] = size(r2); % out-of-sample size and num. of assets

% 1. get covariences and means
sample_covariance = zeros(N,N,T2); % Sample Covariance matrix (sigma hat)
sample_mean = zeros(N,T2); % Sample mean / Expected returns(mu_hat)
for t = 1:T2
    sample_mean(:,t) = mean(r(t:t+T1-1,:));
    sample_covariance(:,:,t) = cov(r(t:t+T1-1,:)); 
end

% general covarience and mean for all T1+T2 values
general_covarience = cov(r); % sigma, size(25,25)
general_mean = mean(r)'; % mu, size(25,1)


% 2. sigma_tilde = (T/(T-N-2))sigma_hat
sigma_tilde = zeros(N,N,T2);
sigma_tilde = T1/(T1-N-2)*sample_covariance;
size(sigma_tilde);


% 3. calculate weight of 1/N and scaled maximum likelyhood
weight_scaled_ML = zeros(N,1,T2); % weight_slash, the weight of scaled maximun likelyhood,size(1,25,T2)
for t=1:T2
    weight_scaled_ML(:,:,t) = 1/gamma * inv(sample_covariance(:,:,t))*sample_mean(:,t);
end

% consider if model is ideal or has transaction cost
transaction_cost = zeros(T2,1); % will be subtracted from each portfolio return. return
if with_transaction_cost
    for t=2:T2
        dw = weight_scaled_ML(:,:,t) - weight_scaled_ML(:,:,t-1); %size(25,1)
        cost = transaction_cost_ratio * abs(dw);
        transaction_cost(t) = sum(cost);
    end
end

weight_1_over_n = ones(N,1); % weight_e, the weight of 1/N. size(1,25)
weight_1_over_n = weight_1_over_n / N;


% 4. calculate pi_1, pi_2, and delta with them
weight_star = inv(general_covarience)*general_mean / gamma; % true optimal portfolio rule, inv(sigma)*mu/gamma, size(1,25)

pi_1 = (weight_1_over_n-weight_star)'*general_covarience*(weight_1_over_n-weight_star);

% Let pi_2 = E(pi_2_list)
pi_2_list = zeros(T2,1);
for t=1:T2
    pi_2_list(t,:) = (weight_scaled_ML(:,:,t)-weight_star)'*general_covarience*(weight_scaled_ML(:,:,t)-weight_star);
end

pi_2 = mean(pi_2_list);
delta = pi_1/(pi_1+pi_2); % combination coefficient


% 5. calculate portfoilo combination 
weight_comb = zeros(N,1,T2); % portfoilo combination, size(25,1,T2)
for t=1:T2
    weight_comb(:,:,t) = (1-delta)*weight_1_over_n + delta*weight_scaled_ML(:,:,t); % portfoilo combination, w_c_hat, size(25,1,T2)
end


% 6. utility
utility = zeros(T2,1); % utility of portfoilo combination(5)
for t=1:T2
    utility(t) = rf2(t) + general_mean'*weight_comb(:,:,t) - gamma/2*weight_comb(:,:,t)'*general_covarience*weight_comb(:,:,t);
end
avg_utility = mean(utility);


% 7. return
portfolio_return = zeros(T2,1);
for t=1:T2
    portfolio_return(t) = rf2(t) + weight_comb(:,:,t)'*r2(t,:)';
end
portfolio_return = portfolio_return - transaction_cost;
avg_return = mean(portfolio_return);


% 8. standard deviation
std_deviation = std(portfolio_return);


% 9. sharpe ratio
sharpe_ratio = mean(portfolio_return) / std_deviation;

end


