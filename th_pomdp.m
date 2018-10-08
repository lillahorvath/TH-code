function pomdp = th_pomdp(pomdp, data_dir, pomdp_dir, plot_pomdp)

% This function creates the PoMDP formulation of the treasure hunt task
% based on its graph representation (nodes,edges) for two treasures. If the 
% plot flag is true, is additionally visualizes the main PoMDP features as 
% shown in Supplementary Fig. 1 of Horvath et al. (201X) XXX.
%
% Inputs
%       pomdp:      structure with fields
%        .d:        grid-world dimensionality
%       data_dir:   group-level data directory
%       pomdp_dir:  name of the pomdp figures directory 
%       plot_pomdp: plot flag
%
% Outputs
%       pomdp:      input structure with additional fields
%        .nn:       number of nodes
%        .N:        set of nodes
%        .E:        set of edges
%        .CN        Cartesian coordinates of nodes
%        .CE        Cartesian coordinates of egdes
%        .A:        set of actions
%        .A_s:      position dependent action set indices
%        .L1:       L1 distance map 
%        .P1:       observation probability map for one target
%        .P2:       observation probability map for two targets
%        .O:        state dependent observation likelihood 
%
% Copyright (C) Dirk Ostwald, Lilla Horvath  
% ------------------------------------------------------------------------- 

% -------------------------------------------------------------------------
%                             initialization 
% -------------------------------------------------------------------------

% unpack input structure
d  = pomdp.d;                                                               
nn = d^2;                                                                   % number of grid-world cells

% -------------------------------------------------------------------------
%                             PoMDP formulation  
% -------------------------------------------------------------------------

% ****************************** graph nodes ******************************

% initialize
N     = NaN(nn,2);                                                          
n_idx = 1;                                                                  

% cycle over the grid-world rows
for i = 1:d
    
    % cycle over the grid-world columns
    for j = 1:d
        
        % set row and column indices
        N(n_idx,1) = i;
        N(n_idx,2) = j;
        
        % update counter
        n_idx = n_idx + 1;
    end
end
   
% ****************************** graph edges ******************************

% initialize
E = zeros(size(N,1));                                                       

% specify square grid-world edges row-wise
for i = 1:(d^2-1)
    
    % right-most column graph nodes
    if mod(i,d) == 0              
        E(i,i+d) = 1;
    
    % last-row column graph nodes
    elseif i >= (d^2 - d)
        E(i,i+1) = 1;
    
    % standard case
    else
        E(i,i+1) = 1;
        E(i,i+d) = 1;
    end
end

% specify symmetric connections
E = E + E';

% *********** conversion to Cartesian coordinate representation ***********

% initialize
CN = NaN(d^2,2);                                                       

% cycle over nodes
for m = 1:size(CN,1)
   
    % convert the matrix graph nodes to Cartesian coords
   CN(m,:) = th_mat2cart(N(m,:),d);
end

% initialize connection array in matrix index space
src_node_mat = NaN(sum(sum(E)),2);
tgt_node_mat = NaN(sum(sum(E)),2);
src_node_crt = NaN(sum(sum(E)),2);
tgt_node_crt = NaN(sum(sum(E)),2);

% cycle over adjacency matrix rows and columns
minidx = 1;
for i = 1:size(E,1)
    for j = 1:size(E,2)
        
        % there exists a connection between M coordinates n and m
        if E(i,j) == 1
            src_node_mat(minidx,:) = N(i,:);
            tgt_node_mat(minidx,:) = N(j,:);
            minidx                 = minidx + 1;    
       end     
    end
end

% cycle over connections
for i = 1:size(src_node_mat,1)
    
    % convert edge coordinates to list of Cartesian coordinates
    src_node_crt(i,:) = th_mat2cart(src_node_mat(i,:),d);
    tgt_node_crt(i,:) = th_mat2cart(tgt_node_mat(i,:),d);
end

% concatenate for output
CE = [src_node_crt tgt_node_crt];

% ******************************** actions ********************************

% create set of grid position dependent actions
GNA = cell(nn,1);                                                           % graph node actions (adjacency list)
A   = [-5 +1 +5 -1];                                                        % action set
A_s = zeros(nn,4);                                                          % position dependent action sets

% cycle over grid-world nodes
for i = 1:numel(GNA)
    
    % identify node indices reachable with one step from node i 
    GNA{i} = find(E(:,i) == 1);
   
    % cycle over actions and determine their agent centered meaning
    for j = 1:length(GNA{i})
        if GNA{i}(j) == i + A(1)
            A_s(i,1) = 1;
        elseif GNA{i}(j) == i + A(2)
            A_s(i,2) = 1;
        elseif GNA{i}(j) == i + A(3)
            A_s(i,3) = 1;
        elseif GNA{i}(j) == i + A(4)
            A_s(i,4) = 1;
        end
    end
end


% ******************** observation likelihood function ********************

% L1 distance map
% -------------------------------------------------------------------------

% ---------------------- l1 distances - single target ---------------------

% initialize L1 distance map
L1 = NaN(nn,nn); 

% cycle over agent positions
for i = 1:nn
    
    % cycle over targets locations - note the reverse row/column order to
    % facilitate row-wise linear indexing 
    for j = 1:nn
        
        % evaluate L1 distance
        L1(i,j) = pdist([N(j,:);N(i,:)], 'cityblock');
    end
end

% ------------------ combined l1 distances - two targets ------------------

% initialize L1 distance map
L1_2 = NaN(nn,nn,nn);

% cycle over agent positions
for i = 1:nn
    
    % cycle over the first target location
    for j = 1:nn
        
        % cycle over the second target location
        for k = 1:nn
            
            % evaluate L1 distance
            L1_2(i,j,k) = min(pdist([N(j,:);N(i,:)], 'cityblock'), pdist([N(k,:);N(i,:)], 'cityblock')) + pdist([N(j,:);N(k,:)], 'cityblock');
        end
    end
end
 
% observation sampling probabilities
% -------------------------------------------------------------------------

% --------- observation sampling probability map - single target ----------

% initialize arrays
beta_0 = NaN(nn,1);                                                         % beta parameter array
beta_1 = NaN(nn,1);                                                         % beta parameter array
P1     = NaN(nn,nn);                                                        % probability map array

% cycle over target locations
for i = 1:nn
    
    % evaluate beta parameters
    beta_0(i) = 1 + (0.5/(max(L1(i,:)) - 1));                               %\beta_0
    beta_1(i) = -(0.5/(max(max(L1(i,:)) - 1)));                             %\beta_1
    
    % evaluate observation sampling probability map
    P1(i,:) = beta_0(i) + beta_1(i).*L1(i,:);                               % row = target location, column = agent position
end

% ----- combined observation sampling probability map - two targets -------

% initialize
P2 = NaN(nn,nn,nn);                                                         % probability map array

% cycle over target 1 locations
for i = 1:nn
    
    % cycle over target 2 locations
    for j = 1:nn
        
        p_1  = reshape(P1(i,:),d,d);                                        % target 1 specific probability maps
        p_2  = reshape(P1(j,:),d,d);                                        % target 2 specific probability maps
        
        p_12 = NaN(d,d);                                                    % initialize combined probability map
        for k = 1:d
            for h = 1:d
                % evaluate combined probability map
                p_12(k,h) = join_p(p_1(k,h),p_2(k,h));                      
            end
        end
        
        % set values 
        P2(i,j,:) = p_12(:);                                                % row = target 1 location, column = target 2 location, page = agent position
    end
end

% state-dependent observation likelihood function p(o|s)
% -------------------------------------------------------------------------

% initialize
O{1} = NaN(nn,nn,nn,4);                                                     % joint observation likelihood 
O{2} = NaN(nn,nn,4);                                                        % marginal observation likelihood 

% -------------- joint observation likelihood - two targets ---------------

% cycle over agent positions
for p_a = 1:nn
    
    % cycle over target 1 locations
    for l_t1 = 1:nn
        
        % cycle over target 2 locations
        for l_t2 = 1:nn
            
             % no target found -> combined probability map 
            O{1}(p_a,l_t1,l_t2, find(A_s(p_a,:))) = P2(l_t1, l_t2, p_a); 

            % check if any action increases the L1 distance from any target
            av_act_idx = find(~isnan(O{1}(p_a,l_t1,l_t2,:)));

            % cycle over available actions
            for a = 1:length(av_act_idx)

                % node index resulting from taking action A(av_act_idx(a))
                act_o = p_a + A(av_act_idx(a));  

                % if action outcome node increases L1 distance to both target 1 and to target 2 set light bar probability to zero
                if L1(act_o,l_t1) >= L1(p_a,l_t1) && L1(act_o,l_t2) >= L1(p_a,l_t2)
                    O{1}(p_a,l_t1,l_t2,av_act_idx(a)) = 0;
                end
            end
        end
    end
end

% -------- marginal observation likelihood function - single target -------

% cycle over agent positions
for p_a = 1:nn
    
    % cycle over target locations  
    for l_t = 1:nn
        
        % one target found -> single target probability map
        O{2}(p_a,l_t,find(A_s(p_a,:))) = P1(l_t,p_a); 

        % check if any action increases the L1 distance from any target
        av_act_idx = find(~isnan(O{2}(p_a,l_t,:)));

        % cycle over available actions
        for a = 1:length(av_act_idx)

            % node index resulting from taking action A(av_act_idx(a))
            act_o = p_a + A(av_act_idx(a));

            % if action outcome node increases L1 distance to target set light bar probability to zero
            if L1(act_o,l_t) >= L1(p_a,l_t) 
                O{2}(p_a,l_t,av_act_idx(a)) = 0;
            end
        end
    end
end

% ************************* set output structure **************************

pomdp.nn  = nn;                                                             % grid-world number of nodes
pomdp.N   = N;                                                              % set of node labels
pomdp.E   = E;                                                              % set of edges
pomdp.CN  = CN;                                                             % Cartesian coordinates of nodes
pomdp.CE  = CE;                                                             % Cartesian coordinates of egdes
pomdp.A   = A;                                                              % action set
pomdp.A_s = A_s;                                                            % position dependent action set indices
pomdp.L1  = L1;                                                             % L1 distance function 
pomdp.P1  = P1;                                                             % observation probability map for one target 
pomdp.P2  = P2;                                                             % observation probability map for two targets
pomdp.O   = O;                                                              % state dependent observation likelihood function

% *********************** visualize PoMDP features ************************

if plot_pomdp
    
    % create directory containing the pomdp figures
    % ---------------------------------------------------------------------
    if exist(fullfile(data_dir, pomdp_dir), 'dir')                          % delete if pre-existent and create a new one
        rmdir(fullfile(data_dir, pomdp_dir),'s')
        mkdir(fullfile(data_dir, pomdp_dir))
    else
        mkdir(fullfile(data_dir, pomdp_dir))                                % create if non-existent
    end
    
    % PoMDP task formulation model components
    % ---------------------------------------------------------------------
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % TO DO: KEEP ONE OF THE FOLLOWING TWO VIS. VERSIONS!!! 
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % ***
    % (1) LARGE: Dirk's print version
    % ***
    
    % ----------------------- figure initialization -----------------------
    sfig1_v1 = figure;
    set( sfig1_v1                ,           ...
         'DefaultFigureColor'    , [1 1 1] , ...
         'DefaultTextInterpreter','Latex'  , ...
         'DefaultAxesFontSize'   , 14      , ...
         'DefaultAxesFontName'   ,'Arial'  , ...
         'DefaultAxesLineWidth'  , 1             )

    mapcol = grey_color_scheme;                                             % color
    mapcol = mapcol/255; 

    % -------------------- subplot 1: nodes and edges ---------------------
    subplot(2,4,1);
    hold on
    plot(CN(:,1), CN(:,2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize',8)

    for i = 1:size(CN,1)
        text(CN(i,1)+0.2, CN(i,2)+0.2, num2str(i), 'FontSize', 10)
    end

    for i = 1:size(CE,1)
        plot([CE(i,1) CE(i,3)],[CE(i,2) CE(i,4)],'k-')
    end
    
    xlim([min(CN(:,1))-.1, max(CN(:,1))+.1])
    ylim([min(CN(:,2))-.1, max(CN(:,2))+.1])
    axis square
    axis off
    title('$N$')

    % -------------------- subplot 2: adjacency matrix --------------------
    ax1 = subplot(2,4,2);
    imagesc(E)
    colormap(ax1,mapcol(2:end,:))
    title('$E$')
    axis square

    % -------------- subplot 3: L1 distances - single target --------------
    ax1 = subplot(2,4,3);
    imagesc(L1)
    colormap(ax1, mapcol(2:20:end,:))
    title('$\ell_1(j,k)$')
    axis square

    % -- subplot 4: observation sampling probability map - single target --
    ax1 = subplot(2,4,4);
    imagesc(P1, [0.5 1])
    colormap(ax1,mapcol(2:end,:))
    title('$p_{j,k}$')
    axis square

    % -- subplot 5: combined observation sampling probability map - two targets --
    ax1 = subplot(2,4,5);
    imagesc(P2(:,:,1), [0.5 1])
    colormap(ax1,mapcol(2:end,:))
    title('$p_{1,k,l}$')
    axis square

    % -- subplot 6: combined observation sampling probability map - two targets --
    ax1 = subplot(2,4,6);
    imagesc(P2(:,:,14), [0.5 1])
    colormap(ax1,mapcol(2:end,:))
    title('$p_{14,k,l}$')
    axis square

    % -- subplot 7: marginal observation likelihood function - single target --
    ax1 = subplot(2,4,7);
    doi = reshape(O{2}(:,22,3),d,d)';                                       % data of interest
    doi(isnan(doi)) = -100;                                                 % replace NaNs with large negative value for visualization purposes
    imagesc(1:d,1:d,doi, [.4 1])
    colormap(ax1,mapcol([1,2:20:end],:))
    axis square
    title('$\pi_{3}\left(j,1,k,22\right)$')

    % ------- subplot 8: likelihood function example - two targets --------
    ax1 = subplot(2,4,8);
    doi = [];
    doi = reshape(O{1}(:,10,22,3),d,d)';                                    % data of interest
    doi(isnan(doi)) = -100;                                                 % replace NaNs with large negative value for visualization purposes
    imagesc(1:d,1:d,doi, [.4 1])
    colormap(ax1,mapcol([1,2:20:end],:))
    axis square
    title('$\pi_{3}\left(j,0,10,22\right)$')

    % --------------------- maximize and print figure ---------------------
    set(sfig1_v1, 'Position', get(0, 'Screensize')); 
    set(sfig1_v1, 'PaperPositionMode', 'auto');
    print('-dbmp', '-r600', fullfile(data_dir, pomdp_dir, 'SFig1_v1'))
    close

    % ***
    % (2) SMALL: Lilla's export_fig version
    % ***
    
    % ----------------------- figure initialization -----------------------
    
    % invoke Yair Altman's export_fig toolbox (usage: save figures with minimal white space)
    addpath(genpath([pwd '\export_fig']))
    
    sfig1_v2 = figure;
    set( sfig1_v2                ,           ...
         'DefaultFigureColor'    , [1 1 1] , ...
         'DefaultTextInterpreter','Latex'  , ...
         'DefaultAxesFontSize'   , 5       , ...
         'DefaultAxesFontName'   ,'Arial'  , ...
         'DefaultAxesLineWidth'  , 0.4           )


    mapcol = grey_color_scheme;                                             % color
    mapcol = mapcol/255; 

    % -------------------- subplot 1: nodes and edges ---------------------
    subplot(2,4,1);
    hold on
    plot(CN(:,1), CN(:,2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize',2)

    for i = 1:size(CN,1)
        text(CN(i,1)+0.2, CN(i,2)+0.2, num2str(i), 'FontSize', 5, 'Interpreter', 'none')
    end

    for i = 1:size(CE,1)
        plot([CE(i,1) CE(i,3)],[CE(i,2) CE(i,4)],'k-', 'Linewidth', 0.3)
    end
    
    xlim([min(CN(:,1))-.1, max(CN(:,1))+.1])
    ylim([min(CN(:,2))-.1, max(CN(:,2))+.1])
    axis square
    axis off
    title('$N$', 'FontSize', 7)

    % -------------------- subplot 2: adjacency matrix --------------------
    ax1 = subplot(2,4,2);
    imagesc(E)
    colormap(ax1,mapcol(2:end,:))
    title('$E$', 'FontSize', 7)
    axis square

    % -------------- subplot 3: L1 distances - single target --------------
    ax1 = subplot(2,4,3);
    imagesc(L1)
    colormap(ax1, mapcol(2:20:end,:))
    title('$\ell_1(j,k)$', 'FontSize', 7)
    axis square

    % -- subplot 4: observation sampling probability map - single target --
    ax1 = subplot(2,4,4);
    imagesc(P1, [0.5 1])
    colormap(ax1,mapcol(2:end,:))
    title('$p_{j,k}$', 'FontSize', 7)
    axis square

    % -- subplot 5: combined observation sampling probability map - two targets --
    ax1 = subplot(2,4,5);
    imagesc(P2(:,:,1), [0.5 1])
    colormap(ax1,mapcol(2:end,:))
    title('$p_{1,k,l}$', 'FontSize', 7)
    axis square

    % -- subplot 6: combined observation sampling probability map - two targets --
    ax1 = subplot(2,4,6);
    imagesc(P2(:,:,14), [0.5 1])
    colormap(ax1,mapcol(2:end,:))
    title('$p_{14,k,l}$', 'FontSize', 7)
    axis square

    % -- subplot 7: marginal observation likelihood function - single target --
    ax1 = subplot(2,4,7);
    doi = reshape(O{2}(:,22,3),d,d)';                                       % data of interest
    doi(isnan(doi)) = -100;                                                 % replace NaNs with large negative value for visualization purposes
    imagesc(1:d,1:d,doi, [.4 1])
    colormap(ax1,mapcol([1,2:20:end],:))
    axis square
    title('$\pi_{3}\left(j,1,k,22\right)$', 'FontSize', 7)

    % ------- subplot 8: likelihood function example - two targets --------
    ax1 = subplot(2,4,8);
    doi = [];
    doi = reshape(O{1}(:,10,22,3),d,d)';                                    % data of interest
    doi(isnan(doi)) = -100;                                                 % replace NaNs with large negative value for visualization purposes
    imagesc(1:d,1:d,doi, [.4 1])
    colormap(ax1,mapcol([1,2:20:end],:))
    axis square
    title('$\pi_{3}\left(j,0,10,22\right)$', 'FontSize', 7)
    
    % ---------------------------- save figure ----------------------------
    sfig1_v2_name = fullfile(data_dir, pomdp_dir, ['SFig1_v2', '.pdf']);
    export_fig(sfig1_v2, sfig1_v2_name);
    close
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % END OF TO DO
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % TO DO: DO WE WANT THE FIG BELOW IN THE PAPER?
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % node dependent available observations & actions
    % ---------------------------------------------------------------------
    
    % ----------------------- figure initialization -----------------------
    h = figure;
    set(h, 'Color', [1 1 1])
    
    % ----------------------------- visualize -----------------------------
    for n = 1:nn
        subplot(d,d,n)
        bar(1:4,A_s(n,:),.7, 'FaceColor', [.8 .8 .8])
        set(gca, 'xtick', 1:4, 'xticklabel', {'$N$', '$E$', '$S$', '$W$'}, 'FontSize', 12, 'TickLabelInterpreter', 'Latex')
        set(gca, 'ytick', 0:1, 'yticklabel', {'$N$', '$A$'}, 'FontSize', 12, 'TickLabelInterpreter', 'Latex')
        title(['$' num2str(n) '$'], 'Interpreter', 'Latex', 'FontSize', 16)
        xlim([0 5])
        ylim([-.2 1.2])
        box off
    end
    
    % --------------------- maximize and print figure ---------------------
    set(h, 'Position', get(0,'Screensize')); 
    set(h, 'PaperPositionMode', 'auto');
    print('-dbmp', '-r600', fullfile(data_dir, pomdp_dir, 'pomdp_act_obs'))
    close
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % END OF TO DO
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
end
end

function p_12 = join_p(p_1,p_2)

% This function evaluates weights w_1, w_2 for two Bernoulli distribution
% parameters p_1, p_2 in [0.5,1] such that the weighted sum of p_1 and p_2
% using w_1, w_2, i.e.
%                      p_12 = w_1*p_2 + w_2*p_2
% falls into the interval [.5,1]. 
%
% Inputs
%       p_1,p_2: Bernoulli distribution parameters in [0.5,1]
%          
% Outputs
%       w_1,w_2: Bernoulli distribution parameter weights in [0,1]
%
% Copyright (C) Dirk Ostwald, Lilla Horvath
% -------------------------------------------------------------------------

if p_1 >= p_2
    p_12 = p_2 + (p_1 - p_2)*p_1;
else
    p_12 = p_1 + (p_2 - p_1)*p_2;
end

end
