function [a_1, a_2, a_3] = th_taskvis_BB_H(taskvis)

% This function visualizes the observable and unobservable data of a
% treasure hunt task instantiation of the hybrid agent models as shown in 
% Supplementary Figs. 8-9 of Horvath et al. (201X) XXX. 
%
% Inputs
%       vistask: structure with fields
%        .pomdp: problem structure
%        .y:     observable data cell array
%        .x:     unobservable data cell array
%        .B:     belief state
%
% Outputs
%       a_1:     figure handle - attempt 1
%       a_2:     figure handle - attempt 2
%       a_3:     figure handle - attempt 3
%
% Copyright (C) Dirk Ostwald, Lilla Horvath
% -------------------------------------------------------------------------

% initialize
% -------------------------------------------------------------------------

% unpack input structure
y  = taskvis.y;                                                             % observable data
x  = taskvis.x;                                                             % unobservable data
B  = taskvis.B;                                                             % belief state
d  = taskvis.pomdp.d;                                                       % grid world dimensionality
CN = taskvis.pomdp.CN;                                                      % Cartesian node coordinates
CE = taskvis.pomdp.CE;                                                      % Cartesian edge coordinates

ne   = numel(y);                                                            % number of episodes
nrow = 4*ne;                                                                % number of plot row                                                         ; % number of plot rows   

% evaluate valid attempts and attempt step limit
sl  = NaN(1,3);
va  = 0;
for a = 1:numel(y)
    sl(a) = size(y{a},1);
    if ~isempty(y{a})
        va = va + 1;
    end
end

% number of plot columns
ncol = max(sl)+1;   

% set default figure values 
fig_pos     = get(0,'defaultfigureposition');
fig_pos_adj = [fig_pos(1:2), 13, 13/fig_pos(3)*fig_pos(4)];
set(0, 'defaultfigurecolor', [1 1 1], 'defaultfigureunits', 'centimeters', 'defaultfigureposition', fig_pos_adj)                                                         

% colors
grey    = grey_color_scheme;
grey_cm = grey./255;
ob_col  = {grey_cm(90,:),grey_cm(210,:)}; 
ag_col  = ice_color_scheme;
ag_col  = ag_col(130,:);
tgt_col = [255, 87, 148 ]./255;

% attempt figure handles
a_1 = [];
a_2 = [];
a_3 = [];

% visualize agent behavior
% -------------------------------------------------------------------------

% cycle over attempts
for a = 1:va
    
    if a == 1
        a_1 = figure;
    elseif a == 2
        a_2 = figure;
    elseif a == 3
        a_3 = figure;
    end
    
    % number of trials to plot (including final trial outcome)
    T = sum(~isnan(y{a}(:,1)));
        
    % cycle over trials 
    for t = 1:T
        
        % information
        p_a = y{a}(t,1);                                                    % agent position
        t_f = y{a}(t,2);                                                    % target found flag
        t_1 = x{a}(t,1);                                                    % target 1 location
        t_2 = x{a}(t,2);                                                    % target 2 location
        
        % --------------------------- grid-world --------------------------
        subplot(nrow,ncol, t+1)
        hold on
        plot(CN(:,1), CN(:,2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize',1)
        for k = 1:size(CE,1)
            plot([CE(k,1) CE(k,3)],[CE(k,2) CE(k,4)],'k-', 'LineWidth', 0.4)
        end

        % ----------------------- agent and targets -----------------------
        
        % no target found and agent not at target location
        if t_f == 0
            plot(CN(p_a,1), CN(p_a,2), 'o', 'MarkerEdgeColor', ag_col, 'MarkerFaceColor', ag_col, 'MarkerSize', 3)
            plot(CN(t_1,1), CN(t_1,2), 'o', 'MarkerEdgeColor', tgt_col, 'MarkerFaceColor', tgt_col, 'MarkerSize', 3)
            plot(CN(t_2,1), CN(t_2,2), 'o', 'MarkerEdgeColor', tgt_col, 'MarkerFaceColor', tgt_col, 'MarkerSize', 3)
 
        % target 1 found 
        elseif t_f == 1 
            plot(CN(p_a,1), CN(p_a,2), 'o', 'MarkerEdgeColor', ag_col, 'MarkerFaceColor', ag_col, 'MarkerSize', 3)
            plot(CN(t_2,1), CN(t_2,2), 'o', 'MarkerEdgeColor', tgt_col, 'MarkerFaceColor', tgt_col, 'MarkerSize', 3)

        % target 2 found 
        elseif t_f == 2 
            plot(CN(p_a,1), CN(p_a,2), 'o', 'MarkerEdgeColor', ag_col, 'MarkerFaceColor', ag_col, 'MarkerSize', 3)
            plot(CN(t_1,1), CN(t_1,2), 'o', 'MarkerEdgeColor', tgt_col, 'MarkerFaceColor', tgt_col, 'MarkerSize', 3)
         
        % both targets found
        elseif t_f == 3
            plot(CN(p_a,1), CN(p_a,2), 'o', 'MarkerEdgeColor', ag_col, 'MarkerFaceColor', ag_col, 'MarkerSize', 3)
        end
               
        xlim([min(CN(:,1)), max(CN(:,1))])
        ylim([min(CN(:,2)), max(CN(:,2))])
        axis square
        axis off
       
        % ------------------------- observations --------------------------
        if t < T
            subplot(nrow,ncol,t+ncol+1)
            for i = 1:4
                % determine a posible observations
                if ~isnan(y{a}(t,i+2))                                      
                    % draw observation bars
                    switch i                                                
                        case 1
                            line([1 2], [2 2],'LineWidth', 3, 'Color', ob_col{y{a}(t,i+2)+1})
                        case 2
                            line([2 2], [1 2],'LineWidth', 3, 'Color', ob_col{y{a}(t,i+2)+1})
                        case 3
                            line([1 2], [1 1],'LineWidth', 3, 'Color', ob_col{y{a}(t,i+2)+1})
                        case 4
                            line([1 1], [2 1],'LineWidth', 3, 'Color', ob_col{y{a}(t,i+2)+1})
                    end
                end
            end
            axis([0.9 2.1 0.9 2.1])
            axis square
            axis off
        end
                
        % ------------- belief state and action probabilities -------------
        
        % belief state
        spb = subplot(nrow,ncol,t+2*ncol);
        imagesc(reshape(B{a}(t,:)',d,d)', [0 1])
        colormap(spb,  [0 0 0; grey_cm(60:end,:)])
        axis square
        axis off

        % normalized BB-R action valences
        if t < T
            ap            = ones(3,3)*-1;
            ap([4 8 6 2]) = x{a}(t,15:18); 
            ax            = subplot(nrow,ncol,t+3*ncol+1);
            imagesc(ap, [-1 0])
            colormap(ax, [0 0 0; grey_cm(60:end,:)])
            axis square
            axis off
        end
        
        % normalized BB-I action valences
        if t < T
            ap            = zeros(3,3);
            ap([4 8 6 2]) = x{a}(t,19:22); 
            ax            = subplot(nrow,ncol,t+4*ncol+1);
            imagesc(ap, [0 1])
            colormap(ax, [0 0 0; grey_cm(60:end,:)])
            axis square
            axis off
        end
        
        % action probabilities
        if t < T
            ap            = zeros(3,3);
            ap([4 8 6 2]) = x{a}(t,3:6);
            ax            = subplot(nrow,ncol,t+5*ncol+1);
            imagesc(ap, [0 1])
            colormap(ax, [0 0 0; grey_cm(60:end,:)])
            axis square
            axis off
        end
        
    end
    
    % final position dependent attempt posterior belief state 
    if ~isempty(B{a})
        spb = subplot(nrow,ncol,t+1+2*ncol);
        imagesc(reshape(B{a}(t+1,:)',d,d)', [0 1])
        colormap(spb, [0 0 0; grey_cm(60:end,:)])
        axis square
        axis off
    end
end




    
    
