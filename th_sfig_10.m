function th_sfig_10(fig_inps)

% This function visualizes basic properties of the tasks as presented to 
% the participants for quality assurance purposes as shown in Supplementary 
% Fig. 10 of Horvath et al. (201X) XXX. The figure consists of the following 
% 3 subplots:
% 
% (1) number of tasks per joint treasure locations 
% (2) number of tasks per shortest path
% (3) number of attempts per step limit noise
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

% invoke Kelly Kearney's plotboxpos toolbox (usage: return position of a plotted region)
addpath(genpath([pwd '\plotboxpos']))

% unpack necessary fields of input structure
data_dir = fig_inps.data_dir;
br_dir   = fig_inps.br_dir;
p_brn    = fig_inps.p_brn;

% load data of interest
load(fullfile(data_dir, br_dir, p_brn)); 

% initialize figure
sfig10      = figure;
set(sfig10, 'Color', [1, 1, 1], 'Units', 'centimeters');
fig_pos     = get(sfig10,'position');
fig_pos_adj = [fig_pos(1:2), 13, 13/fig_pos(3)*fig_pos(4)];
set(sfig10, 'position', fig_pos_adj)
                                                        
% colors
grey_cols = [ 0,0,0;
              100,100,100;
              200,200,200]./255;       
grey_cm   = grey_color_scheme;
grey_cm   = grey_cm./255;
           
% generate figure
% -------------------------------------------------------------------------

% ----------- (1) number of tasks per joint treasure locations  -----------

qa_1_sum = squeeze(sum(qa_1,1));

qa_1_sp = subplot(1, 14, [1 3]);
hold on
qa_1_is = imagesc(qa_1_sum);
set(qa_1_is, 'AlphaData',  ~isnan(qa_1_sum))
xlim([0, pomdp.nn+1])
ylim([0, pomdp.nn+1])
set(gca,'xtick', pomdp.d:pomdp.d:pomdp.nn, 'ytick', pomdp.d:pomdp.d:pomdp.nn, 'FontName', 'Arial', 'FontSize', 5);
set(gca,'Ydir','reverse')
axis square
xlabel('Treasure location', 'FontSize', 6)
ylabel('Treasure location', 'FontSize', 6)
colormap(qa_1_sp, grey_cm(2:end,:));
%colorbar('south', 'Ticks', []);
pos_qa_1 = plotboxpos;

% ---------------- (2) number of tasks per shortest path  -----------------

qa_2_sum = sum(qa_2,1);

qa_2_sp = subplot(1, 14, [5 10]);
hold on
bar(qa_2_sum(2:end), .7, 'FaceColor', grey_cols(2,:), 'EdgeColor', grey_cols(2,:));
set(gca, 'xtick', (0:1:max_opt_dis), 'ytick', 0:20:(max(qa_2_sum)+5), 'FontName', 'Arial', 'FontSize', 5)
ylabel('Number of tasks', 'FontSize', 6)
xlabel('Optimal step size', 'FontSize', 6)
xlim([0, max_opt_dis+1])
ylim([0,(max(qa_2_sum)+1)])
ax                  = gca;
ax.TickDir          = 'out';
ax.XAxis.TickLength = [0,0];
pos_qa_2 = plotboxpos;

% align the size of the subplots
pos_qa_2([2,4]) = pos_qa_1([2,4]);
set(qa_2_sp, 'Position', pos_qa_2);

% ------------- (3) number of attempts per step limit noise  --------------

qa_3_sum = sum(qa_3,1);

qa_3_sp = subplot(1, 14, [12 14]);
hold on
bar(qa_3_sum, .7, 'FaceColor', grey_cols(2,:), 'EdgeColor', grey_cols(2,:));
set(gca, 'xtick', (0:1:sln_set), 'ytick', 0:50:250, 'FontName', 'Arial', 'FontSize', 5)
ylabel('Number of attempts', 'FontSize', 6)
xlabel('Step limit deviation from optimum', 'FontSize', 6)
xlim([0, sln_set+1])
ylim([0, 255])
ax                  = gca;
ax.TickDir          = 'out';
ax.XAxis.TickLength = [0,0];
pos_qa_3 = plotboxpos;

% align the size of the subplots
pos_qa_3([2,4]) = pos_qa_1([2,4]);
set(qa_3_sp, 'Position', pos_qa_3);

% save figure
% -------------------------------------------------------------------------

sfig10_name = fullfile(data_dir, br_dir, ['SFig10', '.pdf']);
export_fig(sfig10, sfig10_name);
close

end

