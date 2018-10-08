function th_fig_3(fig_inps)

% This function visualizes the main participant and agent behavioral 
% results as shown in Fig. 3 of Horvath et al. (201X) XXX. The figure
% consists of the following 6 subplots:
% 
% (1) participant performance
% (2) simulated agent performance - participant task configurations  
% (3) simulated agent performance - standard task configurations 
% (4) participant BIC score
% (5) group cumulative BIC score
% (6) agent model exceedance probability
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
 
% unpack input
data_dir = fig_inps.data_dir;
br_dir   = fig_inps.br_dir;
p_brn    = fig_inps.p_brn;
a_brn    = fig_inps.a_brn;

% load data of interest
load(fullfile(data_dir, br_dir, p_brn));
load(fullfile(data_dir, br_dir, a_brn));

% initialize figure
fig3        = figure;
set(fig3, 'Color', [1, 1, 1], 'Units', 'centimeters');
fig_pos     = get(fig3,'position');
fig_pos_adj = [fig_pos(1:2), 13, 13/fig_pos(3)*fig_pos(4)];
set(fig3, 'position', fig_pos_adj)
                                                        
% colors
grey_cols  = [ 0,0,0;
               100,100,100;
               200,200,200]./255;  
bb_r       = ice_color_scheme;
bb_r       = bb_r([50,90,130,170],:);
bf         = dawn_color_scheme;
bf         = bf([120,160],:);
bb_i       = dusk_color_scheme;
bb_i       = bb_i(150,:);
bb_h       = teal_color_scheme;
bb_h       = bb_h([130, 190],:);
agent_cols = [bf; bb_r; bb_i; bb_h];

% visualization helper variables
x = [0 numel(subs)+1]; 

% generate figure
% -------------------------------------------------------------------------
          
% ----------------------- (1) participant performance ---------------------

bp_1_mean = mean(bp_1);                                                     % mean presented, solvable and solved tasks

subplot(3, 13, [1 6])
hold on
b = bar(bp_1, .7, 'grouped');
for i = 1:3
    set(b(i),'FaceColor', grey_cols(i,:), 'EdgeColor', grey_cols(i,:))
    plot(x, ones(length(x),1).*bp_1_mean(i), 'Color', grey_cols(i,:), 'LineWidth', 0.5, 'LineStyle', '-')
end
set(gca, 'xtick', 1:2:numel(subs), 'ytick', 0:4:16, 'FontName', 'Arial', 'FontSize', 5)
xlabel('Participant', 'FontSize', 6)
ylabel('Number of tasks', 'FontSize', 6)
xlim(x)
ylim([0 16.5])
ax                  = gca;
ax.TickDir          = 'out';
ax.XAxis.TickLength = [0,0];
%lgd = legend('Presented tasks', 'Solvable tasks', 'Solved tasks', 'Location', 'SouthWest')
%lgd.FontSize = 4;

% --- (2) simulated agent performance - participant task configurations ---

sim_pc         = mean(agents_all(:,1:end));                                 % mean agent performance
agents_all_sem = [];                                                        % sem agent performance 
for v = 1:length(mean(agents_all))
    agents_all_sem = [agents_all_sem, nanstd(agents_all(:,v))/sqrt(length(agents_all(~isnan(agents_all(:,v)),v)))];       
end

subplot(3, 13, [8 9])
hold on
for i = 1:numel(agents)
    bar(i, sim_pc(i), .5, 'FaceColor', agent_cols(i,:), 'EdgeColor', agent_cols(i,:))
end
e         = errorbar(1:numel(mean(agents_all(:,1:end))), mean(agents_all(:,1:end)), agents_all_sem(:,1:end), 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 0.5);
e.CapSize = 0.5;
%set(gca, 'xtick', 1:numel(mean(agents_all(:,1:end))), 'XTickLabel',{'BF-0', 'BF-R', 'BB-R-1', 'BB-R-2', 'BB-R-25', 'BB-R-25-2', 'BB-I', 'BB-H-C', 'BB-H-E'}, 'ytick', 0:4:16, 'FontName', 'Arial', 'FontSize', 5)
%xtickangle(60)
set(gca, 'xtick', [], 'ytick', 0:4:16, 'FontName', 'Arial', 'FontSize', 5)
xlabel('Agent', 'Fontsize', 6)
ylabel('Number of solved tasks', 'Fontsize', 6)
xlim([0,numel(agents)+1])
ylim([0,16])
ax         = gca;
ax.TickDir = 'out';
%ax.XAxis.TickLength = [0,0];

% ---- (3) simulated agent performance - standard task configurations -----

avg_perf_loc = squeeze(nanmean(M_sol,2));                                   % average performance over repetitions
avg_perf_abs = nanmean(avg_perf_loc,1);                                     % average performance over repetitions and locations
sem_perf_abs = std(avg_perf_loc,[],1)./sqrt(size(avg_perf_loc,1));          % sem performance over repetitions and locations

subplot(3,13, [11 12])
hold on
for i = 1:numel(agents)
    bar(i, avg_perf_abs(i), .5, 'FaceColor', agent_cols(i,:), 'EdgeColor', agent_cols(i,:));
end
e         = errorbar(1:numel(agents), avg_perf_abs, sem_perf_abs, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 0.5);
e.CapSize = 0.5;
set(gca, 'xtick', [], 'ytick', 0:0.2:1, 'FontName', 'Arial', 'FontSize', 5)
ylabel('% Solved tasks', 'FontSize', 6)
xlabel('Agent', 'FontSize', 6)
xlim([0 numel(agents)+1])
ylim([0 1.02])
ax         = gca;
ax.TickDir = 'out';

% ----------------------- (4) participant BIC score -----------------------

subplot(3, 13, [14 26])
b = bar(bic_all, .25, 'grouped');
for i = 1:numel(agents)
    set(b(i),'FaceColor', agent_cols(i,:), 'EdgeColor', agent_cols(i,:))
end
set(gca, 'xtick', 1:numel(subs), 'ytick', 0:200:600, 'FontName', 'Arial', 'FontSize', 5)
xlabel('Participant', 'Fontsize', 6)
ylabel('BIC score', 'Fontsize', 6)
ylim([0 inf])
xlim(x)
ax                  = gca;
ax.TickDir          = 'out';
ax.XAxis.TickLength = [0,0];
%lgd = legend('BF-0', 'BF-R', 'BB-R-1', 'BB-R-2', 'BB-R-25', 'BB-R-25-2', 'BB-I', 'BB-H-C', 'BB-H-E', 'Location', 'NorthEast', 'Orientation','vertical', 'Linestyle', 'none') 
%lgd.FontSize = 5;
box off

% -------------------- (5) group cumulative BIC score ---------------------

sum_model_BIC = sum(bic_all);                                               % cumulative BIC

subplot(3, 13, [27 28])
hold on
for i = 1:numel(agents)
    bar(i, sum_model_BIC(i), .5, 'FaceColor', agent_cols(i,:), 'EdgeColor', agent_cols(i,:))
end
set(gca, 'xtick',[], 'ytick', 0:2e3:1e4, 'FontName', 'Arial', 'FontSize', 5)
ylabel('Cumulative BIC score', 'Fontsize', 6)
xlabel('Agent', 'Fontsize', 6)
xlim([0 numel(agents)+1])
ylim([0 inf])
ax         = gca;
ax.TickDir = 'out';

% ---------------- (6) agent model exceedance probability -----------------

subplot(3,13, [30 31])
hold on
for i = 1:numel(agents)
    bar(i, av_pp(i), .5, 'FaceColor', agent_cols(i,:), 'EdgeColor', agent_cols(i,:))
end
set(gca, 'xtick', [], 'ytick', 0:0.2:1, 'FontName', 'Arial', 'FontSize', 5)
ylabel('Exceedance probability', 'Fontsize', 6)
xlabel('Agent', 'Fontsize', 6)
xlim([0 numel(agents)+1])
ylim([0 1.02])
ax         = gca;
ax.TickDir = 'out';

% save figure
% -------------------------------------------------------------------------
fig3_name = fullfile(data_dir, br_dir, ['Fig3', '.pdf']);
export_fig(fig3, fig3_name);
close

end

