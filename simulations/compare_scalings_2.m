function compare_scalings_2(target_id)


%% =============  includes...  =============
addpath('utils/')
import_gmp_lib();

%% =============  Load train data  =============
load(['/home/muhsin/workspaces/novel_dmp_generalization-main/simulations/data/comp_scalings_demo.mat'], 'Timed', 'Pd_data', 'dPd_data', 'ddPd_data');

%% =============  Create/Train GMP  =============
n_dof = size(Pd_data, 1); % number of DoFs
n_kernels = 25;

dmp_classic = DMP_classic(n_dof, n_kernels);
classic_mse = dmp_classic.train(Timed, Pd_data,dPd_data, ddPd_data);

dmp_rot = DMP_rot(n_dof, n_kernels);
rot_mse = dmp_rot.train(Timed, Pd_data, dPd_data, ddPd_data);

dmp_bio = DMP_bio(n_dof, n_kernels);
bio_mse = dmp_bio.train(Timed, Pd_data, dPd_data, ddPd_data);

dmp_bio_plus = DMP_bio_plus(n_dof, n_kernels);
bio_plus_mse = dmp_bio_plus.train(Timed, Pd_data, dPd_data, ddPd_data);

gmp = GMP(n_dof, n_kernels, 1.5);
gmp_mse = gmp.train('LS', Timed/Timed(end), Pd_data);
dmp_pp = DMP_pp(gmp);

y0d = Pd_data(:, 1);
gd = Pd_data(:, end);

y0 = y0d;

if target_id == 1
    g = gd + [-0.5; 0.6; 0];
    view_ = [-103.6, 5];
elseif target_id == 2
    g = gd + [-0.4; 0.2; 0];
    view_ = [-103.6, 5];
elseif target_id == 3
    g = gd + [-1.2; 0.6; -0.3];
    view_ = [170.8, 10.33];
elseif target_id == 4
    g = gd + [0.1; -0.3; 0.3];
    view_ = [170.8, 10.33];
elseif target_id == 5
    g = gd + [0.3; 0.2; -0.4];
    view_ = [170.8, 10.33];

else
    error('Invalid target id');
end

Tf = Timed(end);
dt = 0.005;

[Time_cl, P_data_cl, dP_data_cl, ddP_data_cl] = dmp_classic.generate_trajectory_2(y0, g,gd, Tf, dt);
[Time_rot, P_data_rot, dP_data_rot, ddP_data_rot] = dmp_rot.generate_trajectory_2(y0, g,gd, Tf, dt);

[Time_bio, P_data_bio, dP_data_bio, ddP_data_bio] = dmp_bio.generate_trajectory_2(y0, g,gd, Tf, dt);
[Time_bio_plus, P_data_bio_plus, dP_data_bio_plus, ddP_data_bio_plus] = dmp_bio_plus.generate_trajectory_2(y0, g,gd, Tf, dt);

[Time, P_data, dP_data, ddP_data] = dmp_pp.generate_trajectory_2(y0, g,gd, Tf, dt);



%[Time_cl, P_data_cl, dP_data_cl, ddP_data_cl] = dmp_classic.generate_trajectory(y0, g, Tf, dt);
%[Time_rot, P_data_rot, dP_data_rot, ddP_data_rot] = dmp_rot.generate_trajectory(y0, g, Tf, dt);

%[Time_bio, P_data_bio, dP_data_bio, ddP_data_bio] = dmp_bio.generate_trajectory(y0, g, Tf, dt);
%[Time_bio_plus, P_data_bio_plus, dP_data_bio_plus, ddP_data_bio_plus] = dmp_bio_plus.generate_trajectory(y0, g, Tf, dt);

%[Time, P_data, dP_data, ddP_data] = dmp_pp.generate_trajectory(y0, g, Tf, dt);

%% Accumulate the results
dat = {};

dat{1} = struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, ...
                'Color','blue', 'LineStyle','-', 'DisplayName','DMP-PP');

dat{2} = struct('Time',Time_cl, 'Pos',P_data_cl, 'Vel',dP_data_cl, 'Accel',ddP_data_cl, ...
                'Color','magenta', 'LineStyle',':', 'DisplayName','DMP-C');

dat{5} = struct('Time',Time_rot, 'Pos',P_data_rot, 'Vel',dP_data_rot, 'Accel',ddP_data_rot, ...
                'Color','cyan', 'LineStyle','-.', 'DisplayName','DMP-SS');

dat{4} = struct('Time',Time_bio, 'Pos',P_data_bio, 'Vel',dP_data_bio, 'Accel',ddP_data_bio, ...
                 'Color',[0.2 0.8 0.2], 'LineStyle','-', 'DisplayName','DMP-BI');
            
dat{3} = struct('Time',Time_bio_plus, 'Pos',P_data_bio_plus, 'Vel',dP_data_bio_plus, 'Accel',ddP_data_bio_plus, ...
                'Color',[0.85 0.4 0.1], 'LineStyle','-', 'DisplayName','DMP-BIE');

demo = struct('Time',Timed, 'Pos',Pd_data, 'Vel',dPd_data, 'Accel',ddPd_data, ...
                'Color',0.5*[1 1 1], 'LineStyle',':', 'DisplayName','demo');
dat{6} = demo;

% Plot results
%Trajectories
ax = {};
for i=1:3
    figure;
    ax{1} = subplot(3,1,1); hold on; grid on;
    plot(Tf, g(i), 'LineWidth',2, 'Marker','x', 'MarkerSize',12, 'LineStyle','none', 'Color','red', 'HandleVisibility','off');
    for k=1:length(dat)
        if (isempty(dat{k})), continue; end
        plot(dat{k}.Time, dat{k}.Pos(i,:), 'LineWidth',2.0, 'LineStyle',dat{k}.LineStyle, 'Color',dat{k}.Color, 'DisplayName',dat{k}.DisplayName);
    end
    ylabel('pos [$m$]', 'interpreter','latex', 'fontsize',15);
    title(['DoF $' num2str(i) '$'], 'interpreter','latex', 'fontsize',18);
    legend({}, 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;

    ax{2} = subplot(3,1,2); hold on;  grid on;
    for k=1:length(dat)
        if (isempty(dat{k})), continue; end
        plot(dat{k}.Time, dat{k}.Vel(i,:), 'LineWidth',2.0, 'LineStyle',dat{k}.LineStyle, 'Color',dat{k}.Color);
    end
    ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;

    ax{3} = subplot(3,1,3); hold on;  grid on;
    for k=1:length(dat)
        if (isempty(dat{k})), continue; end
        plot(dat{k}.Time, dat{k}.Accel(i,:), 'LineWidth',2.0, 'LineStyle',dat{k}.LineStyle, 'Color',dat{k}.Color);
    end
    ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;

    for j=1:3
        ax{j}.XLim = ax{j}.XLim + 0.01*(ax{j}.XLim(2) - ax{j}.XLim(1)) * [0 1];
        ax{j}.YLim = ax{j}.YLim + 0.02*(ax{j}.YLim(2) - ax{j}.YLim(1)) * [-1 1];
    end
end

% plot path
if n_dof == 3
    plot_ = @(path, varargin) plot3(path(1,:), path(2,:), path(3,:), varargin{:});
elseif n_dof == 2
    plot_ = @(path, varargin) plot(path(1,:), path(2,:), varargin{:}); 
else
    return
end
fig = figure;
fig.Position(3:4) = [686 616];

ax = subplot(2, 1, 1); hold(ax, 'on');
ax.FontSize = 12;
for k=1:length(dat)
    if (isempty(dat{k})), continue; end
    plot_(dat{k}.Pos, 'LineWidth',2.0, 'LineStyle',dat{k}.LineStyle, 'Color',dat{k}.Color, 'DisplayName',dat{k}.DisplayName);
end
plot_(demo.Pos, 'LineWidth',2.0, 'LineStyle',demo.LineStyle, 'Color',demo.Color, 'HandleVisibility','off');
plot_(y0, 'LineWidth',4, 'Marker','o', 'MarkerSize',10, 'LineStyle','none', 'Color','green', 'HandleVisibility','off');
plot_(g, 'LineWidth',4, 'Marker','x', 'MarkerSize',10, 'LineStyle','none', 'Color','red', 'HandleVisibility','off');
plot_(gd, 'LineWidth',4, 'Marker','x', 'MarkerSize',10, 'LineStyle','none', 'Color',[1 0.6 0.6], 'HandleVisibility','off');
legend({}, 'interpreter','latex', 'fontsize',17, 'Position',[0.1052 0.9344 0.7555 0.0506], 'Orientation','horizontal', 'Box','off', 'NumColumns',2);
xlabel('$X$ [$m$]', 'interpreter','latex', 'fontsize',17);
ylabel('$Y$ [$m$]', 'interpreter','latex', 'fontsize',17);
if (n_dof == 3)
    zlabel('$Z$ [$m$]', 'interpreter','latex', 'fontsize',17);
    view(7.7, 11.6);
end
axis tight;
grid on;
hold off;
axis equal;
view(view_);

ax.Position = [0.1300 0.1100 0.7750 0.8150];


addpath('/home/muhsin/workspaces/plots/altmany-export_fig-3.46.0.0');
% 
% %% Figure 1: Plot x, y, z Comparisons
% fig1 = figure;  % Store handle for export
% titles = {'x', 'y', 'z'};  % Titles for the three subplots
% ax = gobjects(3,1);
% colors = lines(length(dat));  % Assign unique colors to each DMP
% 
% for i = 1:3
%     ax(i) = subplot(3,1,i);
%     hold on;
%     grid on;
%     for k = 1:length(dat)
%         if isempty(dat{k}), continue; end
%         plot(dat{k}.Time, dat{k}.Pos(i,:), 'LineWidth',2.0, ...
%             'Color', colors(k,:), 'DisplayName', dat{k}.DisplayName);
%     end
%     % xlabel('Time', 'interpreter','latex', 'fontsize',17);
%     ylabel(titles{i}, 'interpreter','latex', 'fontsize',17);
%     title(titles{i}, 'interpreter','latex', 'fontsize',1);
%     axis tight;
%     hold off;
% end
% 
% % Create a single legend using the first subplot for placement
% % legend(ax(1), 'interpreter','latex', 'fontsize',12, 'Location','bestoutside');
% 
% % Export the first figure as a PDF
% export_fig(fig1, 'XYZ_Comparison.pdf', '-pdf', '-q100', '-nocrop', '-transparent');


%% Figure 1: Plot x, y, z Comparisons
% fig1 = figure;  % Store handle for export
% fig1.Position(3:4) = [600, 500]; 
% titles = {'X', 'Y', 'Z'};  % Use LaTeX-style labels
% titles_2 = {'', '', 'Time steps'};  % Use LaTeX-style labels
% 
% ax = gobjects(3,1);
% 
% % x_limits = [0, 10];  % Set X-axis limits (adjust as needed)
% y_limits = {[-0.3, 0.65], [-2, 2], [-3, 3]}; 
% 
% 
% 
% for i = 1:3
%     ax(i) = subplot(3,1,i);
%     hold on;
%     grid on;
% 
%     plot(Tf, g(i), 'LineWidth',2, 'Marker','x', 'MarkerSize',12, 'LineStyle','none', 'Color','red', 'HandleVisibility','off');
%     for k = 1:length(dat)
%         if isempty(dat{k}), continue; end
%         plot(dat{k}.Time, dat{k}.Pos(i,:), 'LineWidth', 2.0, ...
%             'Color', dat{k}.Color, 'LineStyle', dat{k}.LineStyle, ...
%             'DisplayName', dat{k}.DisplayName);
%     end
%     ylabel(titles{i}, 'interpreter','latex', 'fontsize',17);
%     xlabel(titles_2{i}, 'interpreter','latex', 'fontsize',17);
% 
% 
%     % title(titles{i}, 'interpreter','latex', 'fontsize',17);
% 
%     % xlim(x_limits);  % Apply x-axis limits
%     ylim(y_limits{i});  % Apply y-axis limits%
% 
%     axis tight;
%     hold off;
% end
%% Figure 1: Plot x, y, z Comparisons
fig1 = figure;  % Store handle for export
fig1.Position(3:4) = [500, 490]; 
titles = {'x (m)', 'y (m)', 'z (m)'};  % Use LaTeX-style labels
titles_2 = {'', '', 'Time steps'};  % Use LaTeX-style labels

ax = gobjects(3,1);

% x_limits = [0, 10];  % Set X-axis limits (adjust as needed)
y_limits = {[-0.3, 0.65], [-2, 2], [-3, 3]}; 

% Define the new target (replace with your actual values)
new_target = [0.5, 1.0, 1.5];  % Example new target values for X, Y, Z
desired_distance = [0.23; -0.23; 0.23];
for i = 1:3
    ax(i) = subplot(3,1,i);
    hold on;
    grid on;

        % Plot the demo goal (g(i)) as a solid red dot
    

    % Plot the data from `dat`
    for k = 1:length(dat)
        if isempty(dat{k}), continue; end
        plot(dat{k}.Time, dat{k}.Pos(i,:), 'LineWidth', 2.0, ...
            'Color', dat{k}.Color, 'LineStyle', dat{k}.LineStyle, ...
            'DisplayName', dat{k}.DisplayName);
        demo_endpoint = dat{k}.Pos(i, end);
    end

    plot(Tf, g(i), 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 6, ...
        'LineStyle', 'none', 'Color', 'green', 'MarkerFaceColor', 'green');
    
    % Plot the new target as a solid green dot
    plot(Tf, demo_endpoint, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 6, ...
        'LineStyle', 'none', 'Color', 'red', 'MarkerFaceColor', 'red');


    offset_position = g(i) + desired_distance(i);
    plot(Tf, offset_position, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 6, ...
        'LineStyle', 'none', 'Color', 'blue', 'MarkerFaceColor', 'blue');


    ylabel(titles{i}, 'interpreter','latex', 'fontsize',16);
    xlabel(titles_2{i}, 'interpreter','latex', 'fontsize',16);

    % Apply y-axis limits
    % ylim(y_limits{i});

    set(gca, 'FontSize', 14);

    axis tight;
    hold off;
end

% Add a legend to the plot
% legend('show');


% Create a single legend using the first subplot for placement
% legend(ax(1), 'interpreter','latex', 'fontsize',12, 'Location','bestoutside');

% Export the first figure as a PDF
export_fig(fig1, 'XYZ_Comparison_with_coupling.pdf', '-pdf', '-q100', '-nocrop', '-transparent');




%% Figure 2: Create a Figure with Only a Legend
% fig2 = figure;  % Store handle for export
% hold on;
% grid on;
% colors = lines(length(dat));  % Reuse unique colors for each DMP
% legendEntries = gobjects(1, length(dat));  % Preallocate legend handles
% 
% for k = 1:length(dat)
%     if isempty(dat{k}), continue; end
%     legendEntries(k) = plot(NaN, NaN, 'LineWidth',2.0, ...
%         'Color', colors(k,:), 'DisplayName', dat{k}.DisplayName);
% end
% 
% % Create and format the legend
% legend(legendEntries, 'interpreter','latex', 'fontsize',12, 'Location','best');
% axis off;  % Hide axes so that only the legend is visible
% hold off;
% 
% % Export the legend-only figure as a PDF
% export_fig(fig2, 'XYZ_Legend.pdf', '-pdf', '-q100', '-nocrop', '-transparent');

% %% Figure 2: Create a Legend in a Row
% fig2 = figure;
% fig2.Position(3:4) = [800, 400]; 
% hold on;
% grid on;
% legendEntries = gobjects(1, length(dat));
% 
% for k = 1:length(dat)
%     if isempty(dat{k}), continue; end
%     legendEntries(k) = plot(NaN, NaN, 'LineWidth', 2.0, ...
%         'Color', dat{k}.Color, 'LineStyle', dat{k}.LineStyle, ...
%         'DisplayName', dat{k}.DisplayName);
% end
% 
% % Adjust legend to be in a single row
% legend(legendEntries, 'interpreter','latex', 'fontsize',17, ...
%     'Orientation', 'horizontal', 'Box', 'off', 'NumColumns', length(dat), ... % One row
%     'Position', [0.05 0.9 0.9 0.05]);  % Adjust position to spread across figure
% 
% axis off;
% hold off;
% 
% % Export the legend-only figure as a PDF
% export_fig(fig2, 'XYZ_Legend_classic.pdf', '-pdf', '-q100', '-nocrop', '-transparent');


%% Figure 2: Create a Legend in a Row with Adjustable Box
fig2 = figure;
fig2.Position(3:4) = [1400,80]; 
% fig2.Position(3:4) = [200,80]; 

hold on;
grid on;

% Initialize legend entries for trajectories
legendEntries = gobjects(1, 9);  % +2 for demo goal and new target

% Plot dummy lines for trajectories
for k = 1:length(dat)
    if isempty(dat{k}), continue; end
    legendEntries(k) = plot(NaN, NaN, 'LineWidth', 2.0, ...
        'Color', dat{k}.Color, 'LineStyle', dat{k}.LineStyle, ...
        'DisplayName', dat{k}.DisplayName);
end

% Add demo goal (solid red dot) to legend
legendEntries(end-1) = plot(NaN, NaN, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 8, ...
    'LineStyle', 'none', 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'DisplayName', 'demo goal');
legendEntries(end-2) = plot(NaN, NaN, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 8, ...
    'LineStyle', 'none', 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'DisplayName', 'offset target');
% 
% legendEntries(end-2) = plot(NaN, NaN, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 8, ...
%     'LineStyle', 'none', 'Color', 'black', 'MarkerFaceColor', 'black', ...
%     'DisplayName', 'starting point');

% Add new target goal (solid green dot) to legend
legendEntries(end) = plot(NaN, NaN, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 8, ...
    'LineStyle', 'none', 'Color', 'green', 'MarkerFaceColor', 'green', ...
    'DisplayName', 'new target');

% % Add new target goal (solid green dot) to legend
% legendEntries(end-4) = plot(NaN, NaN, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 8, ...
%     'LineStyle', 'none', 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
%     'DisplayName', 'offset target');

% Create legend with adjustable box
legendHandle = legend(legendEntries, 'interpreter', 'latex', 'fontsize', 17, ...
    'Orientation', 'horizontal', 'Box', 'on', 'NumColumns', length(legendEntries), ...
    'Position', [0.029900850819275 0.305000006258487 0.935152041002324 0.389999987483025]);  % Adjust position and size of the legend box


% Customize the legend box
set(legendHandle, 'EdgeColor', 'black');  % Set border color of the box
set(legendHandle, 'Color', 'white');     % Set background color of the box

axis off;
hold off;

% Export the legend-only figure as a PDF
export_fig(fig2, 'XYZ_Legend_with_coupling.pdf', '-pdf', '-q100', '-transparent');


%% Define n_dof if not already defined
if ~exist('n_dof', 'var')
    n_dof = 3;  % Set to 3 for 3D plotting. Adjust as necessary.
end

%% Define the plotting function based on degrees of freedom
if n_dof == 3
    plot_ = @(path, varargin) plot3(path(1,:), path(2,:), path(3,:), varargin{:});
elseif n_dof == 2
    plot_ = @(path, varargin) plot(path(1,:), path(2,:), varargin{:});
else
    error('n_dof must be 2 or 3');
end

%% Fig 3: 3D Plot Figure
fig3 = figure;
% fig3.Position(3:4) = [686 616];


fig3.Position(3:4) = [600 480];
ax = subplot(2, 1, 1); 
hold(ax, 'on');
ax.FontSize = 14;

% Plot each trajectory from dat
for k = 1:length(dat)
    if isempty(dat{k}), continue; end
    plot_(dat{k}.Pos, 'LineWidth', 2.0, 'LineStyle', dat{k}.LineStyle, ...
          'Color', dat{k}.Color, 'DisplayName', dat{k}.DisplayName);
end

% Plot additional reference trajectories (demo, y0, g, gd)
plot_(demo.Pos, 'LineWidth', 2.0, 'LineStyle', demo.LineStyle, 'Color', demo.Color, 'HandleVisibility','off');
plot_(y0, 'LineWidth', 4, 'Marker', 'o', 'MarkerSize', 6, 'LineStyle', 'none', 'Color', 'black', 'MarkerFaceColor', 'black', 'HandleVisibility','off');
plot_(gd, 'LineWidth', 4, 'Marker', 'o', 'MarkerSize', 6, 'LineStyle', 'none', 'Color', 'red', 'MarkerFaceColor', 'red', 'HandleVisibility','off');
plot_(g, 'LineWidth', 4, 'Marker', 'o', 'MarkerSize', 6, 'LineStyle', 'none', 'Color', 'green', 'MarkerFaceColor', 'green', 'HandleVisibility','off');


% plot_(g + i, 'LineWidth', 4, 'Marker', 'o', 'MarkerSize', 6, ...
%     'LineStyle', 'none', 'Color', 'blue', 'MarkerFaceColor', 'blue', 'HandleVisibility','off');

% Create an (empty) legend - customize as needed
% legend({}, 'interpreter', 'latex', 'fontsize', 17, 'Position', [0.1052 0.9344 0.7555 0.0506], 'Orientation', 'horizontal', 'Box', 'off', 'NumColumns', 2);

xlabel('x (m)', 'interpreter','latex', 'fontsize', 17);
ylabel('y (m)', 'interpreter','latex', 'fontsize', 17);
if n_dof == 3
    zlabel('z (m)', 'interpreter','latex', 'fontsize', 17);
    view(7.7, 11.6);
end

axis tight;
grid on;
hold off;
axis equal;
view(view_);

% Adjust subplot position if needed
ax.Position = [0.1300 0.1100 0.7750 0.8150];

% Export the 3D plot figure as a PDF
export_fig(fig3, '3D_Path_with_coupling.pdf', '-pdf', '-q100', '-transparent');

% 
% dtw_results = struct();
% dtw_results.DMPpp         = dtw(dat{1}.Pos', demo.Pos');        % DMP++ (GMP based)
% dtw_results.DMP_classic   = dtw(dat{2}.Pos', demo.Pos');        % Classic DMP
% dtw_results.DMP_bio_plus  = dtw(dat{3}.Pos', demo.Pos');        % DMP-bio+
% dtw_results.DMP_bio       = dtw(dat{4}.Pos', demo.Pos');        % DMP-bio
% dtw_results.DMP_rot       = dtw(dat{5}.Pos', demo.Pos');        % DMP-rot
% 
% disp('DTW Deviations:');
% disp(dtw_results);


% Compute Euclidean distances
euclidean_results = struct();
euclidean_results.DMPpp = norm(dat{1}.Pos - demo.Pos, 'fro');  % Frobenius norm for Euclidean distance
euclidean_results.DMP_classic = norm(dat{2}.Pos - demo.Pos, 'fro');
euclidean_results.DMP_bio_plus = norm(dat{3}.Pos - demo.Pos, 'fro');
euclidean_results.DMP_bio = norm(dat{4}.Pos - demo.Pos, 'fro');
euclidean_results.DMP_rot = norm(dat{5}.Pos - demo.Pos, 'fro');
disp('Euclidean Distances:');
disp(euclidean_results);


euclidean_results.DMPpp         = sum(vecnorm(dat{1}.Pos' - demo.Pos', 2, 2));
euclidean_results.DMP_classic   = sum(vecnorm(dat{2}.Pos' - demo.Pos', 2, 2));
euclidean_results.DMP_bio_plus  = sum(vecnorm(dat{3}.Pos' - demo.Pos', 2, 2));
euclidean_results.DMP_bio       = sum(vecnorm(dat{4}.Pos' - demo.Pos', 2, 2));
euclidean_results.DMP_rot       = sum(vecnorm(dat{5}.Pos' - demo.Pos', 2, 2));
disp('Euclidean Distances:');
disp(euclidean_results);


% %% --- Append DTW Deviations to a YAML File in Python-compatible Format ---
% % Define the YAML file name
% yaml_filename = 'dtw_deviations_with_coupling.yaml';
% 
% % Define labels and methods
% labels = {'DMP++', 'DMP_classic', 'DMP_bio_plus', 'DMP_bio', 'DMP_rot'};
% methods = {'DTW Deviations'};
% data_values = [dtw_results.DMPpp, dtw_results.DMP_classic, ...
%                dtw_results.DMP_bio_plus, dtw_results.DMP_bio, ...
%                dtw_results.DMP_rot];
% 
% % Open the file in append mode ('a' opens the file for appending)
% fid = fopen(yaml_filename, 'a');
% if fid == -1
%     error('Cannot open file for appending: %s', yaml_filename);
% end
% 
% % Optional: add a newline to separate appended content from existing content
% fprintf(fid, '\n');
% 
% % Write labels section (YAML list format)
% fprintf(fid, 'labels:\n');
% for i = 1:length(labels)
%     fprintf(fid, '  - %s\n', labels{i});
% end
% 
% % Write methods section (YAML list format)
% fprintf(fid, 'methods:\n');
% for i = 1:length(methods)
%     fprintf(fid, '  - %s\n', methods{i});
% end
% 
% % Write data section (YAML nested list format)
% fprintf(fid, 'data:\n');
% fprintf(fid, '  - [');
% fprintf(fid, '%f', data_values(1));  % First value without comma
% for i = 2:length(data_values)
%     fprintf(fid, ', %f', data_values(i));  % Comma-separated values
% end
% fprintf(fid, ']\n');  % Close list
% 
% % Close the file
% fclose(fid);
% 
% disp('YAML file appended successfully in Python-compatibl

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --- Append DTW and Euclidean Deviations to a YAML File in Python-compatible Format ---
% Define the YAML file name
yaml_filename = 'dtw_euclidean_deviations.yaml';

% Define labels and methods
labels = {'DMP++', 'DMP_classic', 'DMP_bio_plus', 'DMP_bio', 'DMP_rot'};
methods = {'DTW Deviations', 'Euclidean Deviations'};

data_values_dtw = [dtw_results.DMPpp, dtw_results.DMP_classic, ...
                   dtw_results.DMP_bio_plus, dtw_results.DMP_bio, ...
                   dtw_results.DMP_rot];
               
% Compute Euclidean distances
data_values_euclidean = sqrt(sum((reference_trajectory - test_trajectory).^2, 2));

data_values_euclidean_avg = mean(data_values_euclidean, 1); % Average over all comparisons

% Open the file in append mode ('a' opens the file for appending)
fid = fopen(yaml_filename, 'a');
if fid == -1
    error('Cannot open file for appending: %s', yaml_filename);
end

% Optional: add a newline to separate appended content from existing content
fprintf(fid, '\n');

% Write labels section (YAML list format)
fprintf(fid, 'labels:\n');
for i = 1:length(labels)
    fprintf(fid, '  - %s\n', labels{i});
end

% Write methods section (YAML list format)
fprintf(fid, 'methods:\n');
for i = 1:length(methods)
    fprintf(fid, '  - %s\n', methods{i});
end

% Write data section (YAML nested list format)
fprintf(fid, 'data:\n');
fprintf(fid, '  - DTW: [');
fprintf(fid, '%f', data_values_dtw(1));
for i = 2:length(data_values_dtw)
    fprintf(fid, ', %f', data_values_dtw(i));
end
fprintf(fid, ']\n');

fprintf(fid, '  - Euclidean: [');
fprintf(fid, '%f', data_values_euclidean_avg(1));
for i = 2:length(data_values_euclidean_avg)
    fprintf(fid, ', %f', data_values_euclidean_avg(i));
end
fprintf(fid, ']\n');

% Close the file
fclose(fid);

disp('YAML file appended successfully with DTW and Euclidean deviations.');





% %% --- Append DTW deviations to a YAML file ---
% % Define the YAML file name
% yaml_filename = 'dtw_deviations_without.yaml';
% 
% % Open the file in append mode ('a' creates the file if it doesn't exist)
% fid = fopen(yaml_filename, 'a');
% if fid == -1
%     error('Cannot open file for appending: %s', yaml_filename);
% end
% 
% % Create a timestamp for this trial
% % timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
% 
% % Append a new YAML document for this trial.
% fprintf(fid, '---\n');
% fprintf(fid, 'trial:\n');
% fprintf(fid, '  target_id: %d\n', target_id);
% % fprintf(fid, '  timestamp: "%s"\n', timestamp);
% fprintf(fid, '  dtw_deviations:\n');
% fprintf(fid, '    DMP++: %f\n', dtw_results.DMPpp);
% fprintf(fid, '    DMP_classic: %f\n', dtw_results.DMP_classic);
% fprintf(fid, '    DMP_bio_plus: %f\n', dtw_results.DMP_bio_plus);
% fprintf(fid, '    DMP_bio: %f\n', dtw_results.DMP_bio);
% fprintf(fid, '    DMP_rot: %f\n', dtw_results.DMP_rot);
% 
% fclose(fid);
% 


% dtw_results = struct();
% dtw_results.DMPpp         = dtw(dat{1}.Pos', demo.Pos');        % DMP++ (GMP based)
% dtw_results.DMP_classic   = dtw(dat{2}.Pos', demo.Pos');        % Classic DMP
% dtw_results.DMP_bio_plus  = dtw(dat{3}.Pos', demo.Pos');        % DMP-bio+
% dtw_results.DMP_bio       = dtw(dat{4}.Pos', demo.Pos');        % DMP-bio
% dtw_results.DMP_rot       = dtw(dat{5}.Pos', demo.Pos');        % DMP-rot
% 
% disp('DTW Deviations:');
% disp(dtw_results);
% 
% %% --- Save DTW deviations to a YAML file ---
% yaml_filename = 'dtw_deviations.yaml';
% fid = fopen(yaml_filename, 'w');
% if fid == -1
%     error('Cannot open file for writing: %s', yaml_filename);
% end
% 
% fprintf(fid, 'dtw_deviations:\n');
% fprintf(fid, '  DMP++: %f\n', dtw_results.DMPpp);
% fprintf(fid, '  DMP_classic: %f\n', dtw_results.DMP_classic);
% fprintf(fid, '  DMP_bio_plus: %f\n', dtw_results.DMP_bio_plus);
% fprintf(fid, '  DMP_bio: %f\n', dtw_results.DMP_bio);
% fprintf(fid, '  DMP_rot: %f\n', dtw_results.DMP_rot);
% 
% fclose(fid);
% disp(['DTW deviations saved to ', yaml_filename]);
% 
% end 


