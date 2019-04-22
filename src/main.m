%test_result = zeros(2,6); % gamma, T1

%gamma = 1; % relative risk aversion coefficient
for with_transaction_cost = 0:1
    for gamma = 1:2:3
        for T1=60:60:360
            fprintf('Test gamma=%d, T1=%d, transaction cost?:%d\n',gamma,T1,with_transaction_cost);
            [avg_utility ,avg_return, std_deviation ,sharpe_ratio, avg_utility_ML ,avg_return_ML, std_deviation_ML ,sharpe_ratio_ML, avg_utility_1_over_n ,avg_return_1_over_n, std_deviation_1_over_n ,sharpe_ratio_1_over_n] = portfolio_assessment(T1,gamma,with_transaction_cost);
            fprintf('For method in 2.2:\n');
            fprintf('avg_utility\t avg_return\t std_deviation\t sharpe_ratio\n');
            fprintf('%.4f\t\t %.4f\t\t %.4f\t\t\t %.4f\t\t\n',avg_utility ,avg_return, std_deviation ,sharpe_ratio);
            fprintf('For maximum likely hood:\n');
            fprintf('avg_utility\t avg_return\t std_deviation\t sharpe_ratio\n');
            fprintf('%.4f\t\t %.4f\t\t %.4f\t\t\t %.4f\t\t\n',avg_utility_ML ,avg_return_ML, std_deviation_ML ,sharpe_ratio_ML);
            fprintf('For 1/N method:\n');
            fprintf('avg_utility\t avg_return\t std_deviation\t sharpe_ratio\n');
            fprintf('%.4f\t\t %.4f\t\t %.4f\t\t\t %.4f\t\t\n',avg_utility_1_over_n ,avg_return_1_over_n, std_deviation_1_over_n ,sharpe_ratio_1_over_n);
            fprintf('-----------------------------------------------------\n');
        end
    end
end


