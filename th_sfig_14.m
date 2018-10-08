function th_sfig_14(fig_inps)

% This function visualizes the results of the model recovery analysis as 
% shown in Supplementary Fig. 14 of Horvath et al. (201X) XXX.
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
sfig14 = figure;
set(sfig14, 'Color', [1, 1, 1], 'Units', 'centimeters');
fig_pos     = get(sfig14,'position');
fig_pos_adj = [fig_pos(1:2), 13, 13/fig_pos(3)*fig_pos(4)];
set(sfig14, 'position', fig_pos_adj)
                                                        
% colors 
bb_r       = ice_color_scheme;
bb_r       = bb_r([50,90,130,170],:);
bf         = dawn_color_scheme;
bf         = bf([120,160],:);
bb_i       = dusk_color_scheme;
bb_i       = bb_i(150,:);
bb_h       = teal_color_scheme;
bb_h       = bb_h([130, 190],:);
agent_cols = [bf; bb_r; bb_i; bb_h];

% generate figure
% -------------------------------------------------------------------------

avg_reco_loc = squeeze(nanmean(M_rec,2));                                   % average agent log likelihoods over repetitions
avg_reco_abs = squeeze(nanmean(avg_reco_loc,1));                            % average agent log likelihoods over repetitions and locations
avg_dn_loc   = squeeze(nanmean(n_dp,2));                                    % average number of data points over repetitions
avg_dn_abs   = squeeze(nanmean(avg_dn_loc,1));                              % average number of data points over repetitions and locations

% compute average BIC scores over repetitions and locations
mrbic_BF_0         = 2.*-avg_reco_abs(:,1);                                 % agent with no free parameters: BF-0
mrbic_BF_R_to_BB_I = 2.*-avg_reco_abs(:,2:7) + log(avg_dn_abs(:,2:7));      % agents with 1 free parameter: BF-R, BB-R-1, BB-R-2, BB-R-25, BB-R-25-2, BB-I
mrbic_BB_H_C       = 2.*-avg_reco_abs(:,8) + 2.*log(avg_dn_abs(:,8));       % agent with 2 free parameters: BB-H-C
mrbic_BB_H_E       = 2.*-avg_reco_abs(:,9) + 3.*log(avg_dn_abs(:,9));       % agent with 3 free parameters: BB-H-E

% concatenate
avg_reco_bic = [mrbic_BF_0, mrbic_BF_R_to_BB_I, mrbic_BB_H_C, mrbic_BB_H_E];

subplot(3, 9, 1:9);
b = bar(avg_reco_bic, .5, 'grouped');
for i = 1:numel(agents)
    set(b(i),'FaceColor', agent_cols(i,:), 'EdgeColor', agent_cols(i,:))
end
set(gca, 'xtick', 1:numel(agents), 'xticklabel', {'BF-0', 'BF-R', 'BB-R-1', 'BB-R-2', 'BB-R-25', 'BB-R-25-2', 'BB-I', 'BB-H-C', 'BB-H-E'}, 'FontName', 'Arial', 'FontSize', 5)
xlabel('Data generating agent', 'FontName', 'Arial', 'Fontsize', 6)
ylabel('BIC score', 'FontName', 'Arial', 'Fontsize', 6)
ylim([0 inf])
xlim([0 9.5])
ax                  = gca;
ax.TickDir          = 'out';
ax.XAxis.TickLength = [0,0];
box off
%lgd = legend('BF-0', 'BF-R', 'BB-R-1', 'BB-R-2', 'BB-R-25', 'BB-R-25-2', 'BB-I', 'BB-H-C', 'BB-H-E', 'Location', 'NorthEast', 'Orientation','vertical', 'Linestyle', 'none');
%lgd.FontSize = 5;

% save figure
% -------------------------------------------------------------------------

sfig14_name = fullfile(data_dir, br_dir, ['SFig14', '.pdf']);
export_fig(sfig14, sfig14_name);
close

end

