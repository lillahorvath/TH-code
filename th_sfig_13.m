function th_sfig_13(fig_inps)

% This function visualizes the results of the agent simulations with
% standard task configurations. The simulation repeats-averaged performance 
% is displayed for each joint treasure location configuration as shown in 
% Supplementary Fig. 13 of Horvath et al. (201X) XXX. 
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
a_brn    = fig_inps.a_brn;

% load data of interest
load(fullfile(data_dir, br_dir, a_brn)); 

% initialize figure
sfig13 = figure;
set(sfig13, 'Color', [1, 1, 1], 'Units', 'centimeters');
fig_pos     = get(sfig13,'position');
fig_pos_adj = [fig_pos(1:2), 13, 13/fig_pos(3)*fig_pos(4)];
set(sfig13, 'position', fig_pos_adj)
                                                        
% colors 
grey_cm = grey_color_scheme;
grey_cm = grey_cm./255;

% generate figure
% -------------------------------------------------------------------------

avg_perf_loc = squeeze(mean(M_sol,2));                                      % average performance over repetitions
map_perf_loc = NaN(pomdp.nn,pomdp.nn,numel(agents));
for a = 1:numel(agents)
    for t = 1:size(tgtloc,1)
        map_perf_loc(tgtloc(t,1), tgtloc(t,2),a) = avg_perf_loc(t,a);
    end
end

a_names = [{'BF-0'}, {'BF-R'}, {'BB-R-1'}, {'BB-R-2'}, {'BB-R-25'}, {'BB-R-25-2'}, {'BB-I'}, {'BB-H-C'}, {'BB-H-E'}];
for a = 1:numel(agents)
    subplot(3,3,a)
    I = imagesc(map_perf_loc(:,:,a), [0 1]);
    set(I, 'alphadata', ~isnan(map_perf_loc(:,:,a)))
    set(gca, 'xtick', 5:10:25, 'ytick', 5:10:25, 'FontName', 'Arial', 'FontSize', 5)
    axis square
    %set(gca,'Ydir','reverse')
    %grid on
    box off
    title([a_names{a}], 'FontName', 'Arial', 'FontSize', 6)
    colormap(grey_cm(2:end-15,:));
    if a == 1
        colorbar('south', 'Ticks', [])
    end
end

% save figure
% -------------------------------------------------------------------------

sfig13_name = fullfile(data_dir, br_dir, ['SFig13', '.pdf']);
export_fig(sfig13, sfig13_name);
close

end

