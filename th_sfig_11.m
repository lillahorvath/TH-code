function th_sfig_11(fig_inps)

% This function visualizes the supplementary participants behavioral results 
% as shown in Supplementary Fig. 11 of Horvath et al. (201X) XXX. The figure 
% consists of the following 4 subplots:
% 
% (1) group performance across runs
% (2) group performance as a function of the number of attempts required 
%     for task solution
% (3) group performance as a function the length of the shortest path
% (4) group performance per joint treasure locations
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
sfig11      = figure;
set(sfig11, 'Color', [1, 1, 1], 'Units', 'centimeters');
fig_pos     = get(sfig11,'position');
fig_pos_adj = [fig_pos(1:2), 13, 13/fig_pos(3)*fig_pos(4)];
set(sfig11, 'position', fig_pos_adj)
                                                        
% colors
grey_cols = [ 0,0,0;
              100,100,100;
              200,200,200]./255;       
grey_cm   = grey_color_scheme;
grey_cm   = grey_cm./255;
           
% generate figure
% -------------------------------------------------------------------------

% ------------------- (1) group performance across runs -------------------

bp_2_mean = nanmean(bp_2);
bp_2_sem  = [];
for v = 1:length(bp_2_mean)
    bp_2_sem = [bp_2_sem, nanstd(bp_2(:,v))/sqrt(length(bp_2(~isnan(bp_2(:,v)),v)))];       
end

bp_2_sp = subplot(2, 5, [1 2]);
hold on
bar(bp_2_mean*100, .7, 'FaceColor', grey_cols(2,:), 'EdgeColor', grey_cols(2,:));
eb         = errorbar(1:numel(bp_2_mean), bp_2_mean*100, bp_2_sem*100, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 0.5);
eb.CapSize = 0.5;
set(gca, 'xtick', 1:max_run, 'ytick', 0:25:100, 'FontName', 'Arial', 'FontSize', 5)
ylabel('% Solved tasks', 'FontSize', 6)
xlabel('Run', 'FontSize', 6)
ylim([0 101])
xlim([0 5])
ax                  = gca;
ax.XAxis.TickLength = [0,0];
ax.TickDir          = 'out';

% ------------------- (2) group performance per attempts ------------------

ap_1_mean = mean(ap_1);
ap_1_std  = std(ap_1);
ap_1_sem  = [];
for v = 1:length(ap_1_std)
    ap_1_sem = [ap_1_sem, nanstd(ap_1(:,v))/sqrt(length(ap_1(~isnan(ap_1(:,v)),v)))];       
end

subplot(2, 5, [4 5]);
hold on
bar(ap_1_mean*100, .7, 'FaceColor', grey_cols(2,:), 'EdgeColor', grey_cols(2,:));
eb         = errorbar(1:numel(ap_1_mean), ap_1_mean*100, ap_1_sem*100, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 0.5);
eb.CapSize = 0.5;
set(gca, 'xtick', 1:max_attempt, 'ytick', 0:10:40, 'FontName', 'Arial', 'FontSize', 5)
ylabel('% Solved tasks', 'FontSize', 6)
xlabel('Attempt', 'FontSize', 6)
xlim([0 max_attempt+1])
ylim([0 41])
ax                  = gca;
ax.XAxis.TickLength = [0,0];
ax.TickDir          = 'out';

% ---------------- (3) group performance per shortest paths ---------------

bp_3a_mean = nanmean(bp_3a,1);
bp_3a_sem  = [];
for v = 1:length(bp_3a_mean)
    bp_3a_sem = [bp_3a_sem, nanstd(bp_3a(:,v))/sqrt(length(bp_3a(~isnan(bp_3a(:,v)),v)))];       
end

subplot(2, 5, [6 7]);
hold on
bar(bp_3a_mean(2:end)*100, .7, 'FaceColor', grey_cols(2,:), 'EdgeColor', grey_cols(2,:));
eb         = errorbar(1:numel(bp_3a_mean)-1, bp_3a_mean(2:end)*100, bp_3a_sem(2:end)*100, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 0.5);
eb.CapSize = 0.5;
set(gca, 'xtick', (1:2:max_opt_dis), 'ytick', 0:25:100, 'FontName', 'Arial', 'FontSize', 5)
ylabel('% Solved tasks', 'FontSize', 6)
xlabel('Shortest path', 'FontSize', 6)
xlim([0, max_opt_dis+1])
ylim([0 101])
ax                  = gca;
ax.XAxis.TickLength = [0,0];
ax.TickDir          = 'out';

% ----------- (4) group performance per joint treasure locations ----------

bp_3b_mean = squeeze(nanmean(bp_3b,1));
for d = 1:pomdp.nn
    for e = d:pomdp.nn
        if isnan(bp_3b_mean(d,e))
            bp_3b_mean(d,e) = -0.5; 
        end
    end
end

bp_3b_sp = subplot(2, 5, [9 10]);
hold on
bp_3b_is = imagesc(bp_3b_mean);
set(bp_3b_is, 'AlphaData',  ~isnan(bp_3b_mean))
xlim([0, pomdp.nn+1])
ylim([0, pomdp.nn+1])
set(gca,'xtick', 5:5:25, 'ytick', 5:5:25, 'FontName', 'Arial', 'FontSize', 5);
xlabel('Treasure location', 'FontSize', 6)
ylabel('Treasure location', 'FontSize', 6)
set(gca,'Ydir','reverse')
axis square
colormap(bp_3b_sp, grey_cm(2:end,:)); 
%cb = colorbar('south', 'Ticks', []);
%set(cb, 'ylim', [0 1])

% save figure
% -------------------------------------------------------------------------

sfig11_name = fullfile(data_dir, br_dir, ['SFig11', '.pdf']);
export_fig(sfig11, sfig11_name);
close

end

