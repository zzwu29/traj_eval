clear;
close all;

path        = matlab.desktop.editor.getActiveFilename;
this_dir    = path(1: end - length(mfilename) - 2);
cd(this_dir);

tests       = dir([this_dir 'result_*']);
tests_count = length(tests);
fprintf('Number of tests: %d\n', length(tests));

rpe_eval = false;

trans_B2prism = [-0.243656,	 -0.012288,	 -0.328095];   % valid with estimation with valid attitude

evaluate_datasets_one_method()

fprintf("-----------------------------------------------------------------------------------------\n")
for i=1:tests_count
    fprintf("%s\t",ATE_POSE{i,1})
    fprintf("%.4f\n",ATE_POSE{i,2})
end
fprintf("-----------------------------------------------------------------------------------------\n")

close all;