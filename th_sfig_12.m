function th_sfig_12(fig_inps)

% This function visualizes the reaction time data as shown in Supplementary 
% Fig. 12 of Horvath et al. (201X) XXX. The figure consists of the following 
% 4 subplots:
% 
% (1) participant median reaction times
% (2) group mean reaction times for each runs
% (3) group mean reaction times for each attempt
% (4) group mean reaction times per shortest paths
%
% Inputs
%       fig_inps:   structure with fields
%        .data_dir: group-level data directory
%        .br_dir:   behavioral results subdirectory: group-level
%        .p_brn:    name of the behavioral results file - participants
%        .a_brn:    name of the behavioral results file - agents
% 
% Outputs
%       none, saves figure directly to disk
% 
% Copyright (C) Lilla Horvath, Dirk Ostwald
% -------------------------------------------------------------------------

% initialize
% -------------------------------------------------------------------------

% unpack necessary fields of input structure
data_dir = fig_inps.data_dir;
br_dir   = fig_inps.br_dir;
p_brn    = fig_inps.p_brn;

% load data of interest
load(fullfile(data_dir, br_dir, p_brn)); 

% initialize figure
sfig12      = figure;
set(sfig12, 'Color', [1, 1, 1], 'Units', 'centimeters');
fig_pos     = get(sfig12,'position');
fig_pos_adj = [fig_pos(1:2), 13, 13/fig_pos(3)*fig_pos(4)];
set(sfig12, 'position', fig_pos_adj)
                                                        
% colors
grey_cols = [ 0,0,0;
              100,100,100;
              200,200,200]./255;       

% visualization helper variables
x = [0 numel(subs)+1]; 
           
% generate figure
% -------------------------------------------------------------------------

% ----------------- (1) participant median reaction times -----------------
brt_1 = brt_1 * 1000;

subplot(2, 10, [1 6]);
hold on
bar(brt_1, .7, 'FaceColor', grey_cols(2,:), 'EdgeColor', grey_cols(2,:));
plot(x, ones(length(x),1).*mean(brt_1), 'Color', grey_cols(1,:), 'LineWidth', 0.5, 'LineStyle', '-')
set(gca, 'xtick', 1:2:numel(subs), 'ytick', 0:200:(max(brt_1)+20), 'FontName', 'Arial', 'FontSize', 5)
xlabel('Participant', 'FontSize', 6)
ylabel('Reaction time (ms)', 'FontSize', 6)
ylim([0 inf])
xlim(x)
ax                  = gca;
ax.TickDir          = 'out';
ax.XAxis.TickLength = [0,0];

% ---------------- (2) group mean reaction times per runs -----------------

brt_2      = brt_2 * 1000;
brt_2_mean = nanmean(brt_2);
brt_2_sem  = [];
for v = 1:length(brt_2_mean)
    brt_2_sem = [brt_2_sem, nanstd(brt_2(:,v))/sqrt(length(brt_2(~isnan(brt_2(:,v)),v)))];       
end

subplot(2, 10, [8 10]);
hold on
bar(brt_2_mean, .7, 'FaceColor', grey_cols(2,:), 'EdgeColor', grey_cols(2,:));
eb = errorbar(1:numel(brt_2_mean), brt_2_mean, brt_2_sem, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 0.5);
eb.CapSize = 0.5;
set(gca, 'xtick', 1:max_run, 'ytick', 0:200:(max(brt_2_mean)+100), 'FontName', 'Arial', 'FontSize', 5)
ylabel('Reaction time (ms)', 'FontSize', 6)
xlabel('Run', 'FontSize', 6)
ylim([0 inf])
xlim([0 5])
ax                  = gca;
ax.TickDir          = 'out';
ax.XAxis.TickLength = [0,0];

% -------------- (3) group mean reaction times per attempts ---------------

art      = art*1000;
art_mean = mean(art);
art_sem  = [];
for v = 1:length(art_mean)
    art_sem = [art_sem, nanstd(art(:,v))/sqrt(length(art(~isnan(art(:,v)),v)))];       
end

subplot(2, 10, [11 13]);
hold on
bar(art_mean, .7, 'FaceColor', grey_cols(2,:), 'EdgeColor', grey_cols(2,:));
eb         = errorbar(1:numel(art_mean), art_mean, art_sem, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 0.5);
eb.CapSize = 0.5;
set(gca, 'xtick', 1:max_attempt, 'ytick', 0:200:1000, 'FontName', 'Arial', 'FontSize', 5)
ylabel('Reaction time (ms)', 'FontSize', 6)
xlabel('Attempt', 'FontSize', 6)
xlim([0 max_attempt+1])
ylim([0 inf])
ax                  = gca;
ax.TickDir          = 'out';
ax.XAxis.TickLength = [0,0];

% ----------- (4) group mean reaction times per shortest paths ------------

brt_3      = brt_3*1000;
brt_3_mean = nanmean(brt_3,1);
brt_3_sem  = [];
for v = 1:length(brt_3_mean)
    brt_3_sem = [brt_3_sem, nanstd(brt_3(:,v))/sqrt(length(brt_3(~isnan(brt_3(:,v)),v)))];       
end

subplot(2, 10, [15 20]);
hold on
bar(brt_3_mean(2:end), .7, 'FaceColor', grey_cols(2,:), 'EdgeColor', grey_cols(2,:));
eb         = errorbar(1:numel(brt_3_mean)-1, brt_3_mean(2:end), brt_3_sem(2:end), 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 0.5);
eb.CapSize = 0.5;
set(gca, 'xtick', (1:2:max_opt_dis), 'ytick', 0:200:(max(brt_3_mean)+100), 'FontName', 'Arial', 'FontSize', 5)
ylabel('Reaction time (ms)', 'FontSize', 6)
xlabel('Shortest path', 'FontSize', 6)
xlim([0, max_opt_dis+1])
ylim([0 inf])
ax                  = gca;
ax.TickDir          = 'out';
ax.XAxis.TickLength = [0,0];

% save figure
% -------------------------------------------------------------------------

sfig12_name = fullfile(data_dir, br_dir, ['SFig12', '.pdf']);
export_fig(sfig12, sfig12_name);
close

end

