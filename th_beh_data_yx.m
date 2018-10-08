function data_yx = th_beh_data_yx(data_yx)

% This function converts the behavioral data that is stored in events long
% format as specified by BIDS to wide format that is more amenable for
% behavioral modelling. In particular, it creates the cell arrays y and x 
% that store the data observable and unobservable to the decision maker,
% respectively.
%
% Inputs
%       data_yx:   structure with fields
%        .sub_dir: subject functional data subdirectory
%        .sub:     subject ID
%        .sub_run: subject run
%        .n_task:  number of tasks in run
%        .a_max:   maximum number of attempts
%
% Outputs
%       data_yx:   input structure with additional fields
%        .sub_dir: subject data directory
%        .sub:     subject ID
%        .sub_run: subject run
%        .n_task:  number of tasks in run
%        .a_max:   maximum number of attempts
%        .y:       cell array containing the observable data
%        .x:       cell array containing the not observable data
%
% Copyright (C) Lilla Horvath, Dirk Ostwald
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%                             initialization 
% -------------------------------------------------------------------------

% unpack input
data_src = data_yx.sub_dir;
subj     = data_yx.sub;
run      = data_yx.sub_run;
n_task   = data_yx.n_task;
a_max    = data_yx.a_max;

% initialize observable and unobservable run-task-attempt arrays
y = cell(1,length(run));
x = cell(1,length(run));

% -------------------------------------------------------------------------
%                convert BIDS events to y and x cell arrays 
% -------------------------------------------------------------------------

% cycle over runs
for r = 1:length(run)

    % ****************** get the BIDS formatted run data ******************
    
    bids_run   = tdfread(fullfile(data_src, [subj, '_task-th_run-0', num2str(r), '_events.tsv']));
    bids_run_c = struct2cell(bids_run);
    
    % convert cell to matrix (convert strings to numbers)
    bids_run_m = zeros(length(bids_run_c{1}),numel(bids_run_c));
    
    % cycle over the events descriptor fields (row-wise entries of the cell array)
    for k = 1:length(bids_run_c)
        % cycle over events
        for l = 1:length(bids_run_c{k})
            if ischar(bids_run_c{k}(l,:))
                if regexp(bids_run_c{k}(l,:), 'n/a')
                    bids_run_m(l,k) = NaN;                                  % switch n/a strings to numeric NaN variables
                else
                    bids_run_m(l,k) = str2num(bids_run_c{k}(l,:));          % convert strings to numbers
                end
            else 
                bids_run_m(l,k) = bids_run_c{k}(l);
            end
        end
    end
    
    % unpack
    ons       = bids_run_m(:,1);                                            % event onset
    dur       = bids_run_m(:,2);                                            % event duration
    ev_type   = bids_run_m(:,3);                                            % event type
    trial_num = bids_run_m(:,4);                                            % trial number
    pos       = bids_run_m(:,5);                                            % agent position
    obs_n     = bids_run_m(:,6);                                            % observation bar - north
    obs_e     = bids_run_m(:,7);                                            % observation bar - east
    obs_s     = bids_run_m(:,8);                                            % observation bar - south
    obs_w     = bids_run_m(:,9);                                            % observation bar - west
    act       = bids_run_m(:,10);                                           % agent action
    tf        = bids_run_m(:,11);                                           % treasure discovery within attempt
    att_num   = bids_run_m(:,12);                                           % attempt number
    tk        = bids_run_m(:,13);                                           % treasure location known from previous attempt
    sli       = bids_run_m(:,14);                                           % attemtp step limit
    task_num  = bids_run_m(:,15);                                           % task number
    tgt_1     = bids_run_m(:,16);                                           % first treasure location
    tgt_2     = bids_run_m(:,17);                                           % second treasure location
    ops       = bids_run_m(:,18);                                           % optimal step size 
    run_data  = [ons, dur, ev_type, trial_num, pos, obs_n, obs_e, obs_s, obs_w, act, tf, att_num, tk, sli, task_num, tgt_1, tgt_2, ops];
    
    % ************************ convert to y and x *************************
    
    % cycle over tasks
    for t = 1:n_task
        
        % get task data
        tioi      = [];
        task_data = [];
        tioi      = find(run_data(:,15) == t);                              % task index of interest
        task_data = run_data(tioi,:);                                       % task data
        
        % cycle over attempts
        for a = 1:a_max
            
            % if attempt exists
            if ismember(a,task_data(:,12))
               
               % get attempt data
               aioi         = [];
               attempt_data = [];
               aioi         = find(task_data(:,12) == a);                   % attempt index of interest
               attempt_data = task_data(aioi, :);                           % attempt data
               ad_y         = NaN(nanmax(attempt_data(:,4)), 11);           % attempt y data
               ad_x         = NaN(nanmax(attempt_data(:,4)), 3);            % attempt x data
               
               % cycle over trials
               for tr = 1:nanmax(attempt_data(:,4))
                   
                   % get trial data
                   trioi      = [];
                   trial_data = [];
                   trioi      = find(attempt_data(:,4) == tr);              % trial index of interest
                   trial_data = attempt_data(trioi,:);                      % trial data
                   
                   % data observable to the decision maker - y cell array
                   % ------------------------------------------------------
                   
                   % the first row of each trial contains all the relevant
                   % observable and unobservable information (this 
                   % corresponds to the agent position presentation event)
                   ad_y(tr,1)   = trial_data(1,5);                          % agent position
                   ad_y(tr,2)   = trial_data(1,11);                         % target discovery flag
                   ad_y(tr,3:6) = trial_data(1,6:9);                        % observation bars (N-E-S-W)
                   ad_y(tr,7)   = trial_data(1,10);                         % actions
                   
                   % compute RT if participant response on trial
                   if ismember(4,trial_data(:,3)) 
                       ad_y(tr,11) = trial_data(logical(trial_data(:,3)==4),1) - trial_data(logical(trial_data(:,3)==3),1);
                   end
                   
                   % constant throughout attempt
                   ad_y(tr,8) = trial_data(1,13);                           % target location known from previous attempt(s)
                   ad_y(tr,9) = trial_data(1,14);                           % step limit
                   
               end
               
               ad_y(:,10) = ones(nanmax(attempt_data(:,4)),1)*ismember(3,attempt_data(:,11)); % task solved in attempt (0: no sol, 1: sol)
               
               % fill in attempt y data
               y{r}{t}{a} = ad_y;
               
               % data not observable to the decision maker - x cell array
               % ----------------------------------------------------------
               ad_x(:,:) = ones(nanmax(attempt_data(:,4)),1)*attempt_data(1,16:18); % target locations, optimal steps
               
              % fill in attempt y data
               x{r}{t}{a} = ad_x;
               
            % if no attempt - keep empty
            else
               y{r}{t}{a} = [];
               x{r}{t}{a} = [];
           end
        end
        
    end
    
    % clear run variables
    clearvars -except subj run data_src y x data_yx n_task a_max
    
end

% ************************* set output structure **************************
data_yx.y = y;
data_yx.x = x;

end
