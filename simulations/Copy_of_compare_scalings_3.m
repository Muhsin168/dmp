% % % Data for the five trials
% % data = [
% %     17.614975, 21.487912, 22.345357, 20.869680, 22.441478;  % Trial 1
% %     13.714779, 34.319106, 31.491924, 32.955572, 32.217275;  % Trial 2
% %     13.714779, 34.319106, 31.491924, 32.955572, 32.217275;  % Trial 3
% %     13.714779, 34.319106, 31.491924, 32.955572, 32.217275;  % Trial 4
% %     13.714779, 34.319106, 31.491924, 32.955572, 32.217275   % Trial 5
% % ];
% % 
% % % Labels for the methods
% % labels = {'DMP++', 'DMP_classic', 'DMP_bio_plus', 'DMP_bio', 'DMP_rot'};
% % 
% % % Create a figure for the box plot
% % figure;
% % 
% % % Generate the box plot
% % boxplot(data', 'Labels', labels);
% % 
% % % Title and axis labels
% % title('Box Plot of DTW Deviations Across Trials');
% % xlabel('Methods');
% % ylabel('DTW Deviations');
% % 
% % % Show the grid
% % grid on;
% % 
% % % Optionally, you can customize the appearance of the box plot
% % 
% 
% 
% 
% % Original Data for the five trials (as you have already provided)
% % Data for the original and improved methods
% % For example:
% % Original and improved data
% % original_data = [
% %     17.614975, 21.487912, 22.345357, 20.869680, 22.441478;  % Trial 1 (original)
% %     13.714779, 34.319106, 31.491924, 32.955572, 32.217275;  % Trial 2 (original)
% %     13.714779, 34.319106, 31.491924, 32.955572, 32.217275;  % Trial 3 (original)
% %     13.714779, 34.319106, 31.491924, 32.955572, 32.217275;  % Trial 4 (original)
% %     13.714779, 34.319106, 31.491924, 32.955572, 32.217275   % Trial 5 (original)
% % ];
% % 
% % improved_data = [
% %     14.614975, 22.487912, 21.345357, 22.869680, 23.441478;  % Trial 1 (improved)
% %     15.714779, 35.319106, 32.491924, 33.955572, 33.217275;  % Trial 2 (improved)
% %     16.714779, 33.319106, 30.491924, 32.955572, 31.217275;  % Trial 3 (improved)
% %     14.714779, 34.319106, 32.491924, 33.955572, 32.217275;  % Trial 4 (improved)
% %     15.714779, 33.319106, 32.491924, 33.955572, 32.217275   % Trial 5 (improved)
% % ];
% % 
% % % Labels for the methods
% % labels = {'DMP++', 'DMP_classic', 'DMP_bio_plus', 'DMP_bio', 'DMP_rot'};
% % 
% % % Create a larger figure
% % figure('Position', [100, 100, 1600, 500]);  % Larger figure size
% % 
% % % Define subplot positions
% % subplot_width = 0.15;  % Width of each subplot
% % subplot_height = 0.7;  % Height of each subplot
% % horizontal_gap = 0.04; % Horizontal gap between subplots
% % left_margin = 0.05;    % Left margin
% % bottom_margin = 0.2;   % Bottom margin
% % 
% % % Create 5 separate subplots to compare original and improved methods
% % for i = 1:5
% %     % Calculate the position for each subplot
% %     left = left_margin + (i-1) * (subplot_width + horizontal_gap);
% %     bottom = bottom_margin;
% % 
% %     % Create a subplot at the calculated position
% %     subplot('Position', [left, bottom, subplot_width, subplot_height]);
% % 
% %     % Combine the original and improved data for each trial (2 sets per trial)
% %     trial_data = [original_data(i, :); improved_data(i, :)];
% % 
% %     % Generate the box plot for the combined data
% %     boxplot(trial_data', 'Labels', {'Original', 'Improved'});
% % 
% %     % Set the title for each subplot (label it with the trial number)
% %     title(['Trial ' num2str(i)]);
% % 
% %     % Set the y-axis label
% %     ylabel('DTW Deviations');
% % 
% %     % Add grid
% %     grid on;
% % end
% % 
% % % Add overall title for the entire figure
% % sgtitle('Comparison of DTW Deviations for Original vs. Improved Methods');
% 
% 
% % Original and improved data
% original_data = [
%     21.561017, 22.136999, 27.086236, 27.321164, 29.223583;  % Trial 1 (new data)
%     11.762319, 12.081718, 12.284313, 14.909446, 14.333489;  % Trial 2 (new data)
%     41.157986, 49.080477, 46.346014, 52.129014, 47.756615;  % Trial 3 (new data)
%     13.714779, 20.933381, 16.314920, 17.376642, 16.990966;  % Trial 4 (new data)
%     17.614975, 27.232579, 23.835659, 22.341656, 23.819448   % Trial 5 (new data)
% ];
% 
% improved_data = [
%     21.561017, 17.149977, 18.417407, 18.369351, 19.295373;  % Trial 1 (new data)
%     11.762319, 9.361524, 9.766024, 10.026638, 10.423294;    % Trial 2 (new data)
%     41.157986, 35.157953, 33.569829, 35.039856, 34.930976;  % Trial 3 (new data)
%     13.714779, 13.355011, 11.564961, 11.679382, 11.751794;   % Trial 4 (new data)
%     17.614975, 17.282373, 15.996808, 15.017940, 15.988451    % Trial 5 (new data)
% ];
% 
% % Skip Trial 1 of improved data
% improved_data(1, :) = [];  % Remove the first row
% 
% % Labels for the methods
% labels = {'DMP++', 'DMP classic', 'DMP bio plus', 'DMP bio', 'DMP rot'};
% 
% 
% 
% figure('Position', [10, 10, 1200, 300]);  % Larger figure size
% 
% % Define subplot positions
% subplot_width = 0.15;  % Width of each subplot
% subplot_height = 0.6;  % Increased height of each subplot to reduce top margin
% horizontal_gap = 0.045; % Horizontal gap between subplots
% left_margin = 0.05;    % Left margin
% bottom_margin = 0.3;  % Reduced bottom margin to reduce top margin
% 
% % Number of trials (after skipping Trial 1 of improved data)
% num_trials = size(original_data, 1);
% 
% % Create subplots to compare original and improved methods
% for i = 1:num_trials
%     % Calculate the position for each subplot
%     left = left_margin + (i-1) * (subplot_width + horizontal_gap);
%     bottom = bottom_margin;
% 
%     % Create a subplot at the calculated position
%     subplot('Position', [left, bottom, subplot_width, subplot_height]);
% 
%     % Combine the original and improved data for each trial
%     if i == 1
%         % For Trial 1, only plot original data (skip improved data)
%         trial_data = original_data(i, :);
%         boxplot(trial_data', 'Labels', {'Original'});
%     else
%         % For other trials, plot both original and improved data
%         trial_data = [original_data(i, :); improved_data(i-1, :)];
%         boxplot(trial_data', 'Labels', {'Original', 'Improved'});
%     end
% 
%     % Set the y-axis label only for the first subplot
%     if i == 1
%         ylabel('DTW Deviations', 'FontSize', 14);  % Increased font size
%     else
%         % Remove y-axis label for other subplots
%         ylabel('');
%     end
% 
%     % Set the x-axis label with increased font size
%     xlabel(labels{i}, 'FontSize', 14, 'Color', 'blue');  % Increased font size
% 
%     % Increase the font size of the tick labels
%     set(gca, 'FontSize', 12);
% 
%     % Add grid
%     grid on;
% end
% % 
% % % Add overall title for the entire figure and adjust its position
% % sgtitle('Comparison of DTW Deviations for Original vs. Improved Methods (Skipping Trial 1 Improved Data)', ...
% %     'FontSize', 14, 'FontWeight', 'bold', 'Position', [0.5, 0.95]);  % Adjust title position


% Data
% Data
% Data and labels (same as before)
original_data = [
    21.561017, 22.136999, 27.086236, 27.321164, 29.223583;  % Trial 1 (new data)
    11.762319, 12.081718, 12.284313, 14.909446, 14.333489;  % Trial 2 (new data)
    41.157986, 49.080477, 46.346014, 52.129014, 47.756615;  % Trial 3 (new data)
    13.714779, 20.933381, 16.314920, 17.376642, 16.990966;  % Trial 4 (new data)
    17.614975, 27.232579, 23.835659, 22.341656, 23.819448   % Trial 5 (new data)
];

improved_data = [
    21.561017, 17.149977, 18.417407, 18.369351, 19.295373;  % Trial 1 (new data)
    11.762319, 9.361524, 9.766024, 10.026638, 10.423294;    % Trial 2 (new data)
    41.157986, 35.157953, 33.569829, 35.039856, 34.930976;  % Trial 3 (new data)
    13.714779, 13.355011, 11.564961, 11.679382, 11.751794;   % Trial 4 (new data)
    17.614975, 17.282373, 15.996808, 15.017940, 15.988451    % Trial 5 (new data)
];

labels = {'DMP++', 'DMP classic', 'DMP bio plus', 'DMP bio', 'DMP rot'};
data = [original_data; improved_data];
groups = [repmat({'Original'}, size(original_data, 1), 1); repmat({'Improved'}, size(improved_data, 1), 1)];

% Create a figure
figure ('Position', [100, 100, 1200, 200]);

% Adjust subplot positions to reduce gaps
num_plots = size(data, 2); % Number of DMP methods
gap = 0.05; % Gap between subplots
margin = 0.09; % Margin around the figure
width = (1 - 2 * margin - (num_plots - 1) * gap) / num_plots; % Width of each subplot

for i = 1:num_plots
    % Calculate position for each subplot
    left = margin + (i - 1) * (width + gap);
    bottom = margin;
    pos = [left, bottom, width, 1 - 2 * bottom]; % [left, bottom, width, height]

    % Create subplot
    axes('Position', pos);
    boxplot(data(:, i), groups);
    title(labels{i});
    ylabel('Performance');
    % xlabel('Data Type');
end

% Adjust layout
% sgtitle('Comparison of Original and Improved Data for Each DMP Method');

% tiledlayout(1, num_plots, 'TileSpacing', 'tight', 'Padding', 'tight');
% 
% for i = 1:num_plots
%     nexttile;
%     boxplot(data(:, i), groups);
%     title(labels{i});
%     ylabel('Performance');
%     xlabel('Data Type');
% end


% Data
original_data = [
    21.561017, 22.136999, 27.086236, 27.321164, 29.223583;  % Trial 1 (new data)
    11.762319, 12.081718, 12.284313, 14.909446, 14.333489;  % Trial 2 (new data)
    41.157986, 49.080477, 46.346014, 52.129014, 47.756615;  % Trial 3 (new data)
    13.714779, 20.933381, 16.314920, 17.376642, 16.990966;  % Trial 4 (new data)
    17.614975, 27.232579, 23.835659, 22.341656, 23.819448   % Trial 5 (new data)
];

improved_data = [
    17.149977, 18.417407, 18.369351, 19.295373;  % Trial 1 (new data, exclude DMP++)
    9.361524, 9.766024, 10.026638, 10.423294;    % Trial 2 (new data, exclude DMP++)
    35.157953, 33.569829, 35.039856, 34.930976;  % Trial 3 (new data, exclude DMP++)
    13.355011, 11.564961, 11.679382, 11.751794;   % Trial 4 (new data, exclude DMP++)
    17.282373, 15.996808, 15.017940, 15.988451    % Trial 5 (new data, exclude DMP++)
];

% Labels for DMP methods
labels = {'DMP++', 'DMP classic', 'DMP bio plus', 'DMP bio', 'DMP rot'};

% Combine original and improved data for each DMP method
% Note: Original data includes all methods, improved data excludes DMP++
data = [original_data; [nan(size(improved_data, 1), 1), improved_data]];

% Create grouping variable for boxplot
% For DMP++, only 'Original' is used; for others, both 'Original' and 'Improved' are used
groups = [repmat({'Original'}, size(original_data, 1), 1); 
          repmat({'Improved'}, size(improved_data, 1), 1)];

% Create a figure with custom size
figure('Position', [100, 100, 1200, 200]); % [left, bottom, width, height]

% Plot boxplots for each DMP method
num_plots = size(data, 2); % Number of DMP methods
gap = 0.03; % Gap between subplots
margin = 0.055; % Margin around the figure
width = (1 - 2 * margin - (num_plots - 1) * gap) / num_plots; % Width of each subplot

for i = 1:num_plots
    % Calculate position for each subplot
    left = margin + (i - 1) * (width + gap);
    bottom = margin;
    pos = [left, bottom, width, 1 - 2 * bottom]; % [left, bottom, width, height]

    % Create subplot
    axes('Position', pos);
    
    % For DMP++, plot only 'Original' data
    if i == 1
        boxplot(data(1:size(original_data, 1), i), groups(1:size(original_data, 1)));
    else
        boxplot(data(:, i), groups);
    end
    
    title(labels{i});
    ylabel('Performance');
    xlabel('Data Type');
end