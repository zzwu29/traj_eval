tic

ATE_POSE    = cell(tests_count, 2);

%% multi-method in one test compare
for n=1:tests_count
    P_h_ate = evaluate_multi_methods(n, tests(n), method_count, methods, trans_B2prism, gt_time_offset, rpe_eval);
    ATE_POSE(n, :) = {tests(n).name, P_h_ate};
end

save('evaluation_result.mat', 'ATE_POSE');

toc