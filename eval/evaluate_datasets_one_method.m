tic

ATE_POSE    = cell(tests_count, 2);

for n=1:tests_count
% parfor n=1:tests_count
    
    P_h_ate = evaluate_one_method(n, tests(n), trans_B2prism, rpe_eval);
            
    ATE_POSE(n, :) = {tests(n).name, P_h_ate};
end


save('evaluation_result.mat', 'ATE_POSE');
toc