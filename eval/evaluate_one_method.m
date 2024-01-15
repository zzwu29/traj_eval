function P_est_ate = evaluate_one_method(test_id, test_fullname, trans_B2prism, rpe_eval)

fig_wid = 550;
% fontset
set(gca, 'fontname','Arial', 'fontsize', 13);

close all;

fighd = [];

%% get the exp name number
exp_name =  test_fullname.name;
exp_path = [test_fullname.folder '/' test_fullname.name '/'];

gndtr_pos_fn     = [exp_path 'leica_pose.csv'];
pose_est_fn      = [exp_path 'predict_odom.csv'];
% trans_B2prism_fn = [exp_path '../trans_B2prism.csv'];


%% Read the gndtr data from logs

% Position groundtr
gndtr_pos_data = csvread(gndtr_pos_fn,  1, 0);

% First sample time used for offsetting all others
t0_ns = gndtr_pos_data(1, 1);

% pos groundtruthdata
t = (gndtr_pos_data(:, 1) - t0_ns)/1e9;
P = gndtr_pos_data(:, 4:6);
Q = gndtr_pos_data(:, [10, 7:9]);

% Delete the duplicate in position groundtruth data
[~, Px_unq_idx] = unique(P(:, 1));
[~, Py_unq_idx] = unique(P(:, 2));
[~, Pz_unq_idx] = unique(P(:, 3));

[~, Qw_unq_idx] = unique(Q(:, 1));
[~, Qx_unq_idx] = unique(Q(:, 1));
[~, Qy_unq_idx] = unique(Q(:, 2));
[~, Qz_unq_idx] = unique(Q(:, 3));

P_unq_idx = union(union(Px_unq_idx, Py_unq_idx), Pz_unq_idx);
Q_unq_idx = union(union(union(Qw_unq_idx, Qx_unq_idx), Qy_unq_idx), Qz_unq_idx);

P = P(P_unq_idx, :);
Q_test = Q(Q_unq_idx, :);
gt_no_att=0;
if  size(Q_test,1)==1 % gt have no att for evaluate
    gt_no_att=1;
end
t = t(P_unq_idx, :);


%% Read the viralslam estimate data from logs
% SLAM estimate
pose_est_data = csvread(pose_est_fn, 1, 0);
t_est = (pose_est_data(:, 1) - t0_ns)/1e9;
P_est =  pose_est_data(:, 4:6);
Q_est = (pose_est_data(:, [10, 7:9]));
V_est =  pose_est_data(:, 11:13);

% Transform from body frame to the prism
% trans_B2prism = csvread(trans_B2prism_fn, 0, 0);

% Compensate the position estimate with the prism displacement
P_est = P_est + quatconv(Q_est, trans_B2prism);


%% Resample the ground truth data by estimate data sample times

% Note affix rs[x] is for resampled by [x]

% Find the interpolated time stamps
[rsest_pos_itp_idx(:, 1), rsest_pos_itp_idx(:, 2)] = combteeth(t_est, t, 0.1);

% Remove the un-associatable samples
rsest_nan_idx = find(isnan(rsest_pos_itp_idx(:, 1)) | isnan(rsest_pos_itp_idx(:, 2)));

t_est_full = t_est;
P_est_full = P_est;
Q_est_full = Q_est;
V_est_full = V_est;

rsest_pos_itp_idx(rsest_nan_idx, :) = [];
t_est(rsest_nan_idx, :)     = [];
P_est(rsest_nan_idx, :)     = [];
Q_est(rsest_nan_idx, :)     = [];
V_est(rsest_nan_idx, :)     = [];

% interpolate the pos gndtr state
P_rsest = vecitp(P, t, t_est, rsest_pos_itp_idx);
Q_rsest=quatitp_slerp(Q, t, t_est, rsest_pos_itp_idx);

% find the optimal alignment
[rot_align_est, trans_align_est] = traj_align(P_rsest, P_est);

% Align the position estimate
P_est      = (rot_align_est*P_est'      + trans_align_est)';
P_est_full = (rot_align_est*P_est_full' + trans_align_est)';

% Align the orientaton estimate
Q_est      = quatmultiply(rotm2quat(rot_align_est), Q_est);
Q_est_full = quatmultiply(rotm2quat(rot_align_est), Q_est_full);

% Align the velocity estimate
V_est      = (rot_align_est*V_est')';
V_est_full = (rot_align_est*V_est_full')';


% Export the leica transform to a yaml file
fileID = fopen([exp_path 'leica_tf.yaml'], 'w');
fprintf(fileID, ['%%YAML:1.0\n'...
                 'T_W_Wleica: !!opencv-matrix\n'...
                 '  rows: 4\n'...
                 '  cols: 4\n'...
                 '  dt: d\n']);
R_W2L   =  rot_align_est';
t_W2L   = -rot_align_est'*trans_align_est;
T_W2L   = [R_W2L, t_W2L; 0 0 0 1];
T_W2L_str = sprintf(['  data: [ %0.9f, %0.9f, %0.9f, %0.9f,\n'...
                     '          %0.9f, %0.9f, %0.9f, %0.9f,\n'...
                     '          %0.9f, %0.9f, %0.9f, %0.9f,\n'...
                     '          %0.9f, %0.9f, %0.9f, %0.9f ]'],...
    T_W2L(1, 1), T_W2L(1, 2), T_W2L(1, 3), T_W2L(1, 4),...
    T_W2L(2, 1), T_W2L(2, 2), T_W2L(2, 3), T_W2L(2, 4),...
    T_W2L(3, 1), T_W2L(3, 2), T_W2L(3, 3), T_W2L(3, 4),...
    T_W2L(4, 1), T_W2L(4, 2), T_W2L(4, 3), T_W2L(4, 4));
fprintf(fileID, T_W2L_str);
fclose(fileID);

% Note: this transform can transform the leica estimate to the slam local
% frame, which can be convenient if you want to record the simulation on
% rviz

mymap_hex=["023EFF", "1AC938", "E8000B","8B2BE2", "FFC400", "00D7FF"]; % bright6 in seaborn
mymap_hex=["90CAF9", "E57373", "C5E1A5", "FFB74D", "F48FB1", "9E86C9", "FFF176", "E6E6E6"]; % https://cdn.elifesciences.org/author-guide/tables-colour.pdf
mymap=zeros(length(mymap_hex),3);
for it=1:length(mymap_hex)
    mymap(it,:)=hex2rgb(mymap_hex(it));
end


%% Calculate the position and rotation errors

%% Calculate the absolute trajectory error of position estimate
P_est_err     = P_rsest - P_est;
P_est_rmse    = rms(P_est_err);
P_est_ate     = norm(P_est_rmse);


%% Print the result
fprintf('Dataset%2d: %s. Err: P_est_ate: %.4f\n',...
          test_id, exp_name(8:end), P_est_ate);



%% Calculate the maximum time
t_max = max([t; t_est]);


% ba_plot_style = {'linestyle', 'none',...
%                   'marker', 'diamond',...
%                   'markerfacecolor', myorange,...
%                   'markeredgecolor', myorange,...
%                   'markersize', 5};


%% Plot the 3D trajectory
% figpos = [1920 0 0 0] + [0, 480, 630, 400];
figpos = [0 0 0 0] + [0, 580, fig_wid, 400];
figure('position', figpos, 'color', 'w', 'paperpositionmode', 'auto');
fighd = [fighd gcf];
hold on;

% % Plot the signal point to let the legend generator use the line symbol
% plot3(P(1:2, 1), P(1:2, 2), P(1:2, 3),  'linewidth', 3);
% plot3(P_est(1:2, 1),  P_est(1:2, 2),  P_est(1:2, 3), 'linewidth', 3);
% 
% % Plot the full trajectory in '.' style to avoid messy gaps
% plot3(P(:, 1), P(:, 2), P(:, 3), '.', 'markersize', 6);
% plot3(P_est_full(:, 1),  P_est_full(:, 2),  P_est_full(:, 3),...
%       '.', 'markersize', 6);

plot3(P(:, 1),  P(:, 2),  P(:, 3), '-o', 'linewidth', 2, 'markersize',  0.01);
plot3(P_est(:, 1),  P_est(:, 2),  P_est(:, 3), '-o', 'linewidth', 2, 'markersize',  0.01);

xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
grid on;
daspect([1 1 1]);
% view([-21 15]);
ax=gca;
ax.ColorOrder=mymap;
tightfig;
set(gca, 'fontsize', 13);
% lg_hd = legend('Leica', 'LOAM (H)', 'LOAM (V)', 'viralslam');
lg_hd = legend('Groundtruth', 'Pos. estimate');

% Save the plot as .fig as well as .png
saveas(gcf, [exp_path exp_name '_traj.fig']);
img = getframe(gcf);
imwrite(img.cdata, [exp_path exp_name '_traj.png']);



%% Plot the time evolution of position
% figpos = [1920 0 0 0] + [0, 0, 630, 400];
figpos = [0 0 0 0] + [fig_wid, 580, fig_wid, 400];
figure('position', figpos, 'color', 'w');
fighd = [fighd gcf];

subplot(3, 1, 1);
hold on;
axgndtr = plot(t,     P(:, 1),    'linewidth', 4);
axest   = plot(t_est, P_est(:, 1),  'linewidth', 2);
uistack(axgndtr, 'top');
uistack(axest, 'top');
ylabel('X (m)');
grid on;
set(gca, 'fontsize', 13);
xlim([0 t_max]);

ax=gca;
ax.ColorOrder=mymap;


subplot(3, 1, 2);
hold on;
axgndtr = plot(t,     P(:, 2),    'linewidth', 4);
axest   = plot(t_est, P_est(:, 2),  'linewidth', 2);
uistack(axgndtr, 'top');
uistack(axest,   'top');
ylabel('Y (m)');
grid on;
set(gca, 'fontsize', 13);
xlim([0 t_max]);

ax=gca;
ax.ColorOrder=mymap;


subplot(3, 1, 3);
hold on;
axgndtr = plot(t,   P(:, 3),    'linewidth', 3);
axest   = plot(t_est,  P_est(:, 3),     'linewidth', 2);
uistack(axgndtr, 'top');
uistack(axest,   'top');
xlabel('Time (s)');
ylabel('Z (m)');
grid on;
set(gca, 'fontsize', 13);
xlim([0 t_max]);

ax=gca;
ax.ColorOrder=mymap;
tightfig(gcf);
saveas(gcf, [exp_path exp_name '_xyzt.fig']);
% saveas(gcf, [exp_path exp_name '_xyzt.pdf']);
img = getframe(gcf);
imwrite(img.cdata, [exp_path exp_name '_xyzt.png']);



%% Plot the time evolution of position estimation error
% figpos = [1920 0 0 0] + [630, 0, 630, 400];
figpos = [0 0 0 0] + [fig_wid * 2, 580, fig_wid, 400];
figure('position', figpos, 'color', 'w');
fighd = [fighd gcf];

subplot(3, 1, 1);
hold on;
plot(t_est,  P_est_err(:, 1),  'linewidth', 2);
ylabel('X Error (m)');
grid on;
set(gca, 'fontsize', 13);
xlim([0 t_max]);

ax=gca;
ax.ColorOrder=mymap;

subplot(3, 1, 2);
hold on;
plot(t_est,  P_est_err(:, 2),   'linewidth', 2);
ylabel('Y Error (m)');
grid on;
set(gca, 'fontsize', 13);
xlim([0 t_max]);

ax=gca;
ax.ColorOrder=mymap;

subplot(3, 1, 3);
hold on;
plot(t_est,  P_est_err(:, 3), 'linewidth', 2);
xlabel('Time (s)');
ylabel('Z Error (m)');
grid on;
set(gca, 'fontsize', 13);
xlim([0 t_max]);

ax=gca;
ax.ColorOrder=mymap;

ax=gca;
ax.ColorOrder=mymap;
tightfig(gcf);
saveas(gcf, [exp_path exp_name '_xyz_err_t.fig']);
img = getframe(gcf);
imwrite(img.cdata, [exp_path exp_name '_xyz_err_t.png']);



%% Plot the combined time evolution of position estimation error
% figpos = [1920 0 0 0] + [630, 480, 630, 200];
figpos = [0 0 0 0] + [0, 100, fig_wid, 400];
figure('position', figpos, 'color', 'w');
fighd = [fighd gcf];

hold on;
plot(t_est, P_est_err(:, 1),  'linewidth', 2);
plot(t_est, P_est_err(:, 2),  'linewidth', 2);
plot(t_est, P_est_err(:, 3),  'linewidth', 2);
xlabel('Time (s)');
ylabel('Error (m)');
grid on;
set(gca, 'fontsize', 13);
xlim([0 t_max]);

lg_hd = legend('Px error', 'Py error', 'Pz error');

ax=gca;
ax.ColorOrder=mymap;
tightfig(gcf);
saveas(gcf, [exp_path exp_name '_xyz_h_err_t.fig']);
img = getframe(gcf);
imwrite(img.cdata, [exp_path exp_name '_xyz_h_err_t.png']);

%% Calculate the relative position error of position estimate

if rpe_eval == true
%% RPE settings
seg_dist = 20;
seg_num = 5;
err_rate = 0.01;

rpe_all=NaN(size(P_rsest,1),seg_num);
rpe_all_perc=NaN(size(P_rsest,1),seg_num);
condition_names = num2cell([1:seg_num]*seg_dist);

fprintf("Wait the rpe pair to match ......")

for i = 1:seg_num
    fprintf(" "+string(i));
    idx_all=get_dist_idx(P_rsest, seg_dist * i, err_rate);
    if gt_no_att==1
        Q_rsest=Q_est;
    end
    rpe_all(1:size(idx_all,1),i)=rpe(idx_all,Q_rsest,Q_est,P_rsest,P_est);
    rpe_all_perc(1:size(idx_all,1),i)=rpe_all(1:size(idx_all,1),i)./(seg_dist * i).*100;
end

fprintf(" ok!\n");

figpos = [0 0 0 0] + [fig_wid, 100, fig_wid, 400];
figure('position', figpos, 'color', 'w');
fighd = [fighd gcf];
hold on;

h = daboxplot(rpe_all,'outliers',0,'xtlabels', condition_names,'fill',0,'colors',mymap);
ylabel('Translation Error (m)');
xlabel('Distance Traveled (m)');
xl = xlim; xlim([xl(1), xl(2)]);     % make more space for the legend

set(gca, 'fontsize', 13);

ax=gca;
ax.ColorOrder=mymap;
tightfig(gcf);
saveas(gcf, [exp_path exp_name '_rpe_err.fig']);
img = getframe(gcf);
imwrite(img.cdata, [exp_path exp_name '_rpe_err.png']);

figpos = [0 0 0 0] + [fig_wid * 2, 100, fig_wid, 400];
figure('position', figpos, 'color', 'w');
fighd = [fighd gcf];
hold on;

h = daboxplot(rpe_all_perc,'outliers',0,'xtlabels', condition_names,'fill',0,'colors',mymap);
ylabel('Translation Error (%)');
xlabel('Distance Traveled (m)');
xl = xlim; xlim([xl(1), xl(2)]);     % make more space for the legend

set(gca, 'fontsize', 13);

ax=gca;
ax.ColorOrder=mymap;
tightfig(gcf);
saveas(gcf, [exp_path exp_name '_rpe_err_perc.fig']);
img = getframe(gcf);
imwrite(img.cdata, [exp_path exp_name '_rpe_err_perc.png']);
end

end