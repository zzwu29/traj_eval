clear;
close all;

path        = matlab.desktop.editor.getActiveFilename;
this_dir    = path(1: end - length(mfilename) - 2);
cd(this_dir);

tests       = dir([this_dir 'result_*']);
tests_count = length(tests);
fprintf('Number of tests: %d\n', length(tests));

%% important! must consistent with your method numbers
% In one test directory, different from evaluate_all.m, the method results
% must be named as predict_odom_method_(i).csv (from 1 to n)
%
methods=["odom1", "odom2", "odom3"]; 
method_count = length(methods);
fprintf('Number of methods: %d\n', length(methods));

%%

% trans_B2prism can be more than a line
% the line num should be the same as method_count

trans_B2prism = [-0.243656,	 -0.012288,	 -0.328095];   % valid with estimation with valid attitude

gt_time_offset = 0.0;   % gt time offset

rpe_eval = true;

evaluate_datasets_multi_methods()

fprintf("-----------------------------------------------------------------------------------------\n")
for i=1:tests_count
    fprintf("%s\t",ATE_POSE{i,1})
    for j=1:method_count
        fprintf("%.4f\t",ATE_POSE{i,2}(j))
    end
    fprintf("\n")
end
fprintf("-----------------------------------------------------------------------------------------\n")

% close all;