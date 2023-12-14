function P_est_ate_all_method = evaluate_multi_methods(test_id, test_fullname, num_methods,lgd_methods, trans_B2prism, rpe_eval)

fig_wid = 550;

%% RPE settings
seg_dist = 20;
seg_num = 5;
err_rate = 0.01; 
err_rate = 0.02;

% fontset
set(gca, 'fontname','Arial', 'fontsize', 13);

close all;

P_est_ate_all_method=[];

fighd = [];

%% get the exp name number
exp_name =  test_fullname.name;
exp_path = [test_fullname.folder '/' test_fullname.name '/'];

gndtr_pos_fn     = [exp_path 'leica_pose.csv'];
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

%% Plot the 3D trajectory
figpos = [0 0 0 0] + [0, 580, fig_wid, 400];
figure('position', figpos, 'color', 'w', 'paperpositionmode', 'auto');
fighd = [fighd gcf];
hold on;
  
plot3(P(:, 1), P(:, 2), P(:, 3), '-','linewidth', 2);

%% Plot the time evolution of position
figpos = [0 0 0 0] + [fig_wid, 580, fig_wid, 400];
figure('position', figpos, 'color', 'w');
fighd = [fighd gcf];

subplot(3, 1, 1);
hold on;
axgndtr = plot(t,     P(:, 1), 'linewidth', 2);

subplot(3, 1, 2);
hold on;
axgndtr = plot(t,     P(:, 2), 'linewidth', 2);

subplot(3, 1, 3);
hold on;
axgndtr = plot(t,   P(:, 3), 'linewidth', 2);


%% Plot the time evolution of position estimation error
figpos = [0 0 0 0] + [fig_wid * 2, 580, fig_wid, 400];
figure('position', figpos, 'color', 'w');
fighd = [fighd gcf];

subplot(3, 1, 1);
hold on;

subplot(3, 1, 2);
hold on;

subplot(3, 1, 3);
hold on;

%% Calculate the relative position error of position estimate

rpe_all=[];
rpe_all_perc=[];
group_inx=[];

condition_names = num2cell([1:seg_num]*seg_dist);
group_names = mat2cell(lgd_methods,1, ones(size(lgd_methods)));

for m=1:num_methods
    %% Read the viralslam estimate data from logs
    % SLAM estimate
    clear rsest_pos_itp_idx;
    pose_est_fn      = [exp_path 'predict_odom_method_'  char(string(m))  '.csv'];

    pose_est_data = csvread(pose_est_fn, 1, 0);
    t_est = (pose_est_data(:, 1) - t0_ns)/1e9;
    P_est =  pose_est_data(:, 4:6);
    Q_est = (pose_est_data(:, [10, 7:9]));
    V_est =  pose_est_data(:, 11:13);

    % Transform from body frame to the prism
%     trans_B2prism = csvread(trans_B2prism_fn, 0, 0);

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

    rsest_pos_itp_idx(rsest_nan_idx, :) = [];
    t_est(rsest_nan_idx, :)     = [];
    P_est(rsest_nan_idx, :)     = [];
    Q_est(rsest_nan_idx, :)     = [];

    % interpolate the pos gndtr state
    P_rsest = vecitp(P, t, t_est, rsest_pos_itp_idx);

    % find the optimal alignment
    [rot_align_est, trans_align_est] = traj_align(P_rsest, P_est);

    % Align the position estimate
    P_est      = (rot_align_est*P_est'      + trans_align_est)';
    P_est_full = (rot_align_est*P_est_full' + trans_align_est)';

    % Align the orientaton estimate
    Q_est      = quatmultiply(rotm2quat(rot_align_est), Q_est);
    Q_est_full = quatmultiply(rotm2quat(rot_align_est), Q_est);


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


    %% Calculate the position and rotation errors


    %% Calculate the absolute trajectory error of position estimate
    P_est_err     = P_rsest - P_est;
    P_est_rmse    = rms(P_est_err);
    P_est_ate     = norm(P_est_rmse);


    %% Print the result
    fprintf('Test%2d: %s. Method%2d: %s. Err: P_est_ate: %.4f\n',...
              test_id, exp_name(8:end), m, lgd_methods(m), P_est_ate);

    P_est_ate_all_method=[P_est_ate_all_method,P_est_ate];

    %% Calculate the maximum time
    t_max = max([t; t_est]);

    %% Plot the 3D trajectory
    figure(1)
    hold on;

    plot3(P_est(:, 1),  P_est(:, 2),  P_est(:, 3), '-o', 'linewidth', 2, 'markersize',  0.01);


    %% Plot the time evolution of position
    figure(2)
    hold on;

    subplot(3, 1, 1);
    hold on;
    axest   = plot(t_est, P_est(:, 1), '-o', 'linewidth', 2, 'markersize',  0.01);

    subplot(3, 1, 2);
    hold on;
    axest   = plot(t_est, P_est(:, 2), '-o', 'linewidth', 2, 'markersize',  0.01);

    subplot(3, 1, 3);
    hold on;
    axest   = plot(t_est,  P_est(:, 3), '-o', 'linewidth', 2, 'markersize', 0.01);



    %% Plot the time evolution of position estimation error
    figure(3)
    hold on;

    subplot(3, 1, 1);
    hold on;
    plot(t_est,  P_est_err(:, 1), '-o', 'linewidth', 2, 'markersize',  0.01);

    subplot(3, 1, 2);
    hold on;
    plot(t_est,  P_est_err(:, 2), '-o', 'linewidth', 2, 'markersize',  0.01);

    subplot(3, 1, 3);
    hold on;
    plot(t_est,  P_est_err(:, 3), '-o', 'linewidth', 2, 'markersize',  0.01);
    
    if rpe_eval == true
    fprintf("Wait the rpe pair in method %2d to match ......", m);

    rpe_iter_all=NaN(size(P_rsest,1),seg_num);
    rpe_iter_all_perc=NaN(size(P_rsest,1),seg_num);
    
    for i = 1:seg_num
        fprintf(" "+string(i));
        idx_all=get_dist_idx(P_rsest, seg_dist * i, err_rate);
        if gt_no_att==1
            Q_rsest=Q_est;
        end
        rpe_iter=rpe(idx_all,Q_rsest,Q_est,P_rsest,P_est);
        rpe_iter_all(1:size(idx_all,1),i)=rpe(idx_all,Q_rsest,Q_est,P_rsest,P_est);
        rpe_iter_all_perc(1:size(idx_all,1),i)=rpe_iter_all(1:size(idx_all,1),i)./(seg_dist * i).*100;
    end
    
    rpe_all=[rpe_all;rpe_iter_all];
    rpe_all_perc=[rpe_all_perc;rpe_iter_all_perc];
    group_inx=[group_inx  m*ones(1,size(rpe_iter_all,1))];
    
    fprintf(" ok!\n");
    end

end
% ba_plot_style = {'linestyle', 'none',...
%                   'marker', 'diamond',...
%                   'markerfacecolor', myorange,...
%                   'markeredgecolor', myorange,...
%                   'markersize', 5};

mymap_hex=["023EFF", "1AC938", "E8000B","8B2BE2", "FFC400", "00D7FF"]; % bright6 in seaborn
mymap_hex=["90CAF9", "E57373", "C5E1A5", "FFB74D", "F48FB1", "9E86C9", "FFF176", "E6E6E6"]; % https://cdn.elifesciences.org/author-guide/tables-colour.pdf
mymap=zeros(length(mymap_hex),3);
for it=1:length(mymap_hex)
    mymap(it,:)=hex2rgb(mymap_hex(it));
end

lgd_names=lgd_methods;
%% Plot the 3D trajectory
figure(1);
hold on;

xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
grid on;
daspect([1 1 1]);
% view([-21 15]);
set(gca, 'fontsize', 13);
% lg_hd = legend('Leica', 'LOAM (H)', 'LOAM (V)', 'viralslam');
lg_hd = legend(["Groundtruth",lgd_names]);
% set(lg_hd,'box','off')
ax=gca;
ax.ColorOrder=mymap;
tightfig(gcf);
% Save the plot as .fig as well as .png
saveas(gcf, [exp_path exp_name '_traj.fig']);
img = getframe(gcf);
imwrite(img.cdata, [exp_path exp_name '_traj.png']);



%% Plot the time evolution of position
figure(2);
hold on;

subplot(3, 1, 1);
hold on;
ylabel('X (m)');
grid on;
 set(gca, 'fontsize', 13);
xlim([0 t_max]);
ax=gca;
ax.ColorOrder=mymap;

subplot(3, 1, 2);
hold on;
ylabel('Y (m)');
grid on;
 set(gca, 'fontsize', 13);
xlim([0 t_max]);
ax=gca;
ax.ColorOrder=mymap;

subplot(3, 1, 3);
hold on;
xlabel('Time (s)');
ylabel('Z (m)');
grid on;
 set(gca, 'fontsize', 13);
xlim([0 t_max]);
ax=gca;
ax.ColorOrder=mymap;


lg_hd = legend(["Groundtruth",lgd_names]);
tightfig(gcf);
saveas(gcf, [exp_path exp_name '_xyzt.fig']);
% saveas(gcf, [exp_path exp_name '_xyzt.pdf']);
img = getframe(gcf);
imwrite(img.cdata, [exp_path exp_name '_xyzt.png']);



%% Plot the time evolution of position estimation error
figure(3);
hold on;

subplot(3, 1, 1);
hold on;
ylabel('X Error (m)');
grid on;
 set(gca, 'fontsize', 13);
xlim([0 t_max]);
ax=gca;
ax.ColorOrder=mymap;

subplot(3, 1, 2);
hold on;
ylabel('Y Error (m)');
grid on;
 set(gca, 'fontsize', 13);
xlim([0 t_max]);
ax=gca;
ax.ColorOrder=mymap;

subplot(3, 1, 3);
hold on;
xlabel('Time (s)');
ylabel('Z Error (m)');
grid on;
set(gca, 'fontsize', 13);
xlim([0 t_max]);
ax=gca;
ax.ColorOrder=mymap;


lg_hd = legend(lgd_names);
tightfig(gcf);
saveas(gcf, [exp_path exp_name '_xyz_err_t.fig']);
img = getframe(gcf);
imwrite(img.cdata, [exp_path exp_name '_xyz_err_t.png']);

if rpe_eval == true
    
figpos = [0 0 0 0] + [fig_wid * 0.5, 100, fig_wid, 400];
figure('position', figpos, 'color', 'w');
fighd = [fighd gcf];
hold on;

h = daboxplot(rpe_all,'groups',group_inx,'outliers',0,'xtlabels', condition_names,'fill',0,'legend',group_names,'colors',mymap);
ylabel('Translation Error (m)');
xlabel('Distance Traveled (m)');
xl = xlim; xlim([xl(1), xl(2)]); 

set(gca, 'fontsize', 13);

ax=gca;
ax.ColorOrder=mymap;
tightfig(gcf);
saveas(gcf, [exp_path exp_name '_rpe_err.fig']);
img = getframe(gcf);
imwrite(img.cdata, [exp_path exp_name '_rpe_err.png']);

figpos = [0 0 0 0] + [fig_wid * 1.5, 100, fig_wid, 400];
figure('position', figpos, 'color', 'w');
fighd = [fighd gcf];
hold on;

h = daboxplot(rpe_all_perc,'groups',group_inx,'outliers',0,'xtlabels', condition_names,'fill',0,'legend',group_names,'colors',mymap);
ylabel('Translation Error (%)');
xlabel('Distance Traveled (m)');
xl = xlim; xlim([xl(1), xl(2)]);  

set(gca, 'fontsize', 13);

ax=gca;
ax.ColorOrder=mymap;
tightfig(gcf);
saveas(gcf, [exp_path exp_name '_rpe_err_perc.fig']);
img = getframe(gcf);
imwrite(img.cdata, [exp_path exp_name '_rpe_err_perc.png']);

end
end
