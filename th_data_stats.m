function [stats] = th_data_stats(stats)

% This function evaluates basic descriptive statistics of the treasure hunt
% behavioural data of a single participant
%
% Inputs
%       stats:         input structure with fields
%        .sub_run:     subject run
%        .pomdp:       POMDP formulation  
%        .n_task:      number of tasks in run
%        .a_max:       max num of attempts
%        .sp_max:      maximum length of shortest path (=optimal step size)
%        .y:           observable data
%        .x:           unobservable data
%
% Outputs
%       stats:         input structure with additional fields
%        .tl:          treasure location statistics
%        .dist:        distance statistics
%        .run:         run statistics
%        .attempt:     attempt statistics
%        .dn:          decision statistics
%        .ntask:       number of tasks presented
%        .solvable:    number of solvable tasks
%        .solved:      number of solved tasks
%        .tasks:       remaining statistics of interest
%        .step_lim:    sampled step limit noise statistics
%        .RT:          basic RT data statistics
%        .run_rt:      run RT stats
%        .dist_rt:     distance RT stats
%        .attempt_rt: attempt RT stats
%
% Copyright (C) Dirk Ostwald, Lilla Horvath
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%                             initialization 
% -------------------------------------------------------------------------

% unpack input structure
pomdp  = stats.pomdp;                                                       % PoMDP task formulation
run    = stats.sub_run;                                                     % subject run 
n_task = stats.n_task;                                                      % number of tasks in run
a_max  = stats.a_max;                                                       % maximum number of attempts
sp_max = stats.sp_max;                                                      % maximum length of shortest path
y      = stats.y;                                                           % observable data
x      = stats.x;                                                           % unobservable data
nruns  = numel(run);                                                        % number of runs

% experiment task index
tskidx = 1;

% initialize task data array
tasks = cell(1,1);

% -------------------------------------------------------------------------
%                       compute summary statistics 
% -------------------------------------------------------------------------

% *************************** get key task data ***************************

% cycle over runs
for r = 1:nruns
    
    % cycle over tasks in run
    for t = 1:numel(y{r})
        
        % initialize array 
        keydata = [];                                                      
        
        % cycle over attempts
        for a = 1:numel(y{r}{t})
            
            if ~isempty(x{r}{t}{a})
                              
                % treasure (=target) locations, smaller node index first
                tgtloc = sort(x{r}{t}{a}(1,1:2))';
                              
                % optimal number of steps
                optstep = x{r}{t}{a}(1,3);
                
                % available number of steps (=step limit)
                avstep = y{r}{t}{a}(1,9);
                
                % target 1 or target 2 found
                if ~isnan(y{r}{t}{a}(end,2))
                    tf_1_2 = y{r}{t}{a}(end,2); 
                else
                    tf_1_2 = y{r}{t}{a}(end-1,2);
                end
                
                % task solved
                solved = y{r}{t}{a}(1,10);
                    
                % concatenate
                keydata = [keydata [tgtloc; optstep; avstep; tf_1_2 ; solved]];                              
                   
            end
        end
        
        % collect keydata
        tasks{tskidx} = keydata;
             
        % increase task index
        tskidx = tskidx + 1;
       
    end
end

% ************* evaluate basic performance related statistics *************

% task statistics
% -------------------------------------------------------------------------

% initialize 
ntasks     = numel(tasks);                                                  % number of tasks
task_stats = NaN(7,ntasks);                                                 % task statistics array
step_lim_m = NaN(a_max,ntasks);                                             % step limit noise counter array 

% cycle over tasks
for t = 1:ntasks
    task_stats(1,t) = tasks{t}(1,1);                                        % target 1 location
    task_stats(2,t) = tasks{t}(2,1);                                        % target 2 location
    task_stats(3,t) = tasks{t}(3,1);                                        % optimal number of steps
    task_stats(4,t) = max(tasks{t}(4,:));                                   % max available number of steps
    task_stats(5,t) = any(tasks{t}(3,:) <= tasks{t}(4,:));                  % task solvability
    task_stats(6,t) = sum(tasks{t}(end,:),2);                               % task solved
    task_stats(7,t) = size(tasks{t},2);                                     % number of attempts
    
    % step limit stats
    step_lim_m(1:task_stats(7,t),t) = tasks{t}(4,:)-tasks{t}(3,:);
    if t == ntasks
        step_lim_sum = [sum(nansum(step_lim_m==-2)), sum(nansum(step_lim_m==-1)), sum(nansum(step_lim_m==0)), sum(nansum(step_lim_m==1)), sum(nansum(step_lim_m==2))];
    end
end

% attempt statistics
% -------------------------------------------------------------------------
attempt_perf = NaN(a_max,2);                                                % attempt performance array
dat_solvab   = task_stats(5,:);                                             % data array of interest - solvability
dat_solved   = task_stats(6,:);                                             % data array of interest - solved

% cycle over attempts
for i = 1:a_max
    ind               = logical(task_stats(7,:) == i);                      % data of interest index
    attempt_perf(i,1) = sum(dat_solved(ind))/sum(dat_solvab);               % percentage of solved in attempt out of all solvable
    attempt_perf(i,2) = sum(dat_solved(ind))/sum(dat_solved);               % percentage of solved in attempt out of all solved
end

% run statistics
% -------------------------------------------------------------------------
run_stats = NaN(nruns,4);                                                   % run statistics array
j         = 0;                                                              % counter

% cycle over runs
for i = 1:nruns                                                 
    run_stats(i,1) = numel(y{i});                                           % number of tasks in run
    run_stats(i,2) = sum(task_stats(5,j+1:j+numel(y{i})));                  % solvable tasks in run
    run_stats(i,3) = sum(task_stats(6,j+1:j+numel(y{i})));                  % solved tasks in run
    run_stats(i,4) = run_stats(i,3)/run_stats(i,2);                         % percentage of solved tasks (out of solvable)
    
    % update counter
    j = j + numel(y{i});
end

% distance statistics
% -------------------------------------------------------------------------
dist_stats      = NaN(sp_max+1, 5);                                         % distance statistics
dist_stats(:,1) = (0:sp_max)';                                              % possible l1 distances as given by the grid layout

% cycle over the set of optimal step sizes
for i = 1:size(dist_stats,1)
    dist_stats(i,2) = sum(task_stats(3,:) == dist_stats(i,1));                      % number of tasks per opt step size
    dist_stats(i,3) = sum(task_stats(5, find(task_stats(3,:) == dist_stats(i,1)))); % number of solvable tasks per opt step size
    dist_stats(i,4) = sum(task_stats(6, find(task_stats(3,:) == dist_stats(i,1)))); % number of solved tasks per opt step size
    dist_stats(i,5) = dist_stats(i,4)/dist_stats(i,3);                              % percentage of solved tasks (out of solvable)
end

% treasure location statistics 
% -------------------------------------------------------------------------
tl_stats = NaN(pomdp.nn,pomdp.nn,4);                                        % tgt location statistics

% cycle over grid nodes
for i = 1:pomdp.nn
    
    % cycle over grid nodes
    for j = 1:i
        tl              = ismember(task_stats(1:2,:)', sort([i,j]), 'rows'); 
        tl_stats(i,j,1) = sum(tl);                                          % number of tasks per joint treasure locations
        solvable        = task_stats(5,:)';
        tl_stats(i,j,2) = sum(solvable(tl));                                % number of solvable tasks per joint treasure locations
        solved          = task_stats(6,:)';
        tl_stats(i,j,3) = sum(solved(tl));                                  % number of solved tasks per joint treasure locations
        tl_stats(i,j,4) = sum(solved(tl))/sum(solvable(tl));                % percentage of solved tasks (out of solvable)
    end
end

% decision statistics
% -------------------------------------------------------------------------
dec_nums    = NaN(6,1);                                                     % decision stats array
num_a       = 0;                                                            % initialize counter: number of actions
missing_a   = 0;                                                            % initialize counter: number of missing actions
num_a_db    = 0;                                                            % initialize counter: number of actions towards dark bars
num_a_lb    = 0;                                                            % initialize counter: number of actions towards light bars
num_a_db_lb = 0;                                                            % initialize counter: number of actions towards dark bars in the presence of light bars
num_a_lb_db = 0;                                                            % initialize counter: number of actions towards light bars in the presence of dark bars

% cycle over runs
for i = 1:nruns
   
    % cycle over tasks
    for j = 1:n_task
       
        % cycle over attempts
        for k = 1:a_max
            
            % if attempt was carried out
            if ~isempty(y{i}{j}{k})
                dat       = [];
                dat       = y{i}{j}{k};
                dat_o     = dat(:,3:6);
                dat_o     = dat_o(~isnan(dat(:,7)),:);
                dat_a     = dat(:,7);
                dat_a     = dat_a(~isnan(dat_a)); 
                missing_a = missing_a + ((length(dat(:,7))-1) - length(dat_a));
                
                % cycle over trials
                for l = 1:length(dat_a)
                    
                    a_db    = 0;
                    a_lb    = 0;
                    a_db_lb = 0;
                    a_lb_db = 0;
                   
                    if dat_a(l) == 1
                        
                        if dat_o(l,2) == 0
                            a_db = 1;
                            
                            if ismember(1, dat_o(l,:))
                                a_db_lb = 1;
                            end
                        end

                    elseif dat_a(l) == -1
                        
                        if dat_o(l,4) == 0
                            a_db = 1;
                            
                            if ismember(1, dat_o(l,:))
                                a_db_lb = 1;
                            end
                        end
                        
                    elseif dat_a(l) == 5
                        
                        if dat_o(l,3) == 0
                            a_db = 1;
                            
                            if ismember(1, dat_o(l,:))
                                a_db_lb = 1;    
                            end
                        end
                        
                    elseif dat_a(l) == -5
                        
                        if dat_o(l,1) == 0
                            a_db = 1;
                            
                            if ismember(1, dat_o(l,:))
                                a_db_lb = 1;   
                            end
                        end
                        
                    end
                    
                    if a_db == 0
                        a_lb = 1;
                        
                        if ismember(0, dat_o(l,:))
                            a_lb_db = 1;
                        end
                    end
                    
                    num_a       = num_a + numel(a_db);
                    num_a_db    = num_a_db + a_db;
                    num_a_lb    = num_a_lb + a_lb;
                    num_a_db_lb = num_a_db_lb + a_db_lb;
                    num_a_lb_db = num_a_lb_db + a_lb_db;
                    
                end 
            end
        end
    end
end
dec_nums(1) = num_a;
dec_nums(2) = missing_a;
dec_nums(3) = num_a_db;
dec_nums(4) = num_a_lb;
dec_nums(5) = num_a_db_lb;
dec_nums(6) = num_a_lb_db;

% **************** evaluate basic reaction time statistics ****************

% save RT data in numerical array 
% -------------------------------------------------------------------------
nt_max  = sp_max + 2;                                                       % maximum number of trials per attempt: max distance + 2 (step limits)
na_max  = ntasks * a_max;                                                   % maximumal number of attempts per participant
rt_data = NaN(nt_max,na_max);
idx     = 1;

% cycle over runs
for r = 1:numel(y)
  
    % cycle over tasks
    for t = 1:numel(y{r})
   
        % cycle over attempts
         for a = 1:numel(y{r}{t})
          
             % extract reaction time data
             if ~isempty(y{r}{t}{a})
                rt_data(1:size(y{r}{t}{a}(:,11),1),idx) = y{r}{t}{a}(:,11);
                idx = idx + 1;
             end
         end
    end
end

% basic RT statistics
% -------------------------------------------------------------------------
valid_rt_idx = ~isnan(rt_data);                                             % obtain indices of valid reaction times
RT.total     = sum(valid_rt_idx(:));                                        % total number of RTs
RT.trial     = sum(valid_rt_idx,2);                                         % number of RTs as function of trial on attempt
RT.avgtrial  = nanmean(rt_data, 2);                                         % average RT as a function of trial on attempt
RT.stdtrial  = nanstd(rt_data, [], 2);                                      % RT standard deviation as a function of trial on attempt
RT.all       = rt_data(valid_rt_idx);                                       % all RTs
RT.avg       = median(RT.all,1);                                            % median RT

% run RTs
% -------------------------------------------------------------------------
nart = NaN(nruns,n_task);                                                   % number of attempts per runs and tasks

% cycle over runs
for i = 1:nruns 
    
    % cycle over tasks
    for j = 1:n_task
       acrt = [];
       
       % cycle over attempts
       for k = 1:a_max
           acrt = [acrt, ~isempty(y{i}{j}{k})];
       end
       nart(i,j) = sum(acrt);
    end
end

nar     = sum(nart,2);                                                      % number of attempts per run 
run_rts = NaN(nruns, 1);                                                    % aggregate RT per run
j       = 0;                                                                % counter

% cycle over runs
for i = 1:length(nar)
    rt_r       = rt_data(:,j+1:j+nar(i));                                   % all trial RT in run
    v_rt_r     = ~isnan(rt_r);                                              % index: valid rt in run
    run_rts(i) = median(rt_r(v_rt_r));                                      % median of all valid RTs in run
    j          = j + nar(i);                                                % update counter
end

% attempt RTs
% -------------------------------------------------------------------------
attempt_rts = NaN(a_max,1);                                                 % aggregate RT per attempt
rt_a        = NaN(nt_max*na_max, a_max);                                    % all trial RT in attempt
t_to_a      = NaN(sum(sum(nart)),1);                                        % task number to attempts
t           = 1;                                                            % counter
zzz         = 1;                                                            % counter

% cycle over runs
for i = 1:nruns
    z = 1; % counter
    
    % cycle over tasks
    for j = 1:n_task
        zz = 1; % counter
        
        % cycle over attempts
        for k = z:sum(nart(i,1:j)) 
            ind = [];
            dat = [];
            ind = find(isnan(rt_a(:,zz)), 1, 'first');                      % index to write data in rt_a
            dat = rt_data(~isnan(rt_data(:,k)),k);                          % attempt data
            
            rt_a(ind:(ind+length(dat)-1),zz) = dat;
            
            z           = z+1;                                              % update counter 
            zz          = zz+1;                                             % update counter
            t_to_a(zzz) = t;                                                % write task number to attempt
            zzz         = zzz+1;                                            % update counter
        end
        t = t+1; % update counter
    end
end

% cycle over attempts
for i = 1:a_max
    attempt_rts(i) = nanmedian(rt_a(:,i));                                  % get median RT per attempt
end

% distance RTs
% -------------------------------------------------------------------------
rt_d     = cell(sp_max+1,1);                                                % all RTs per shortest path 
dist_rts = NaN(numel(rt_d),1);                                              % aggregate RT per shortest path 

% cycle over tasks
for i = 1:ntasks   
    rt_task                 = [];                                           % RTs of the task attempts of interest
    rt_task                 = rt_data(:,t_to_a == i);
    rt_d{task_stats(3,i)+1} = [rt_d{task_stats(3,i)+1}; rt_task(~isnan(rt_task))];
end

% cycle over optimal step sizes
for i = 1:length(dist_rts) 
    dist_rts(i) = median(rt_d{i});                                          % get median RT per shortest path
end
   
% ************************* set output structure **************************
stats.tl         = tl_stats;                                                % treasure location statistics
stats.dist       = dist_stats;                                              % distance statistics
stats.run        = run_stats;                                               % run statistics
stats.attempt    = attempt_perf;                                            % attempt statistics
stats.dn         = dec_nums;                                                % decision statistics
stats.ntask      = ntasks;                                                  % number of tasks presented 
stats.solvable   = sum(task_stats(5,:));                                    % number of solvable task
stats.solved     = sum(task_stats(6,:));                                    % number of solved tasks
stats.tasks      = tasks;                                                   % remaining statistics of interest
stats.step_lim   = step_lim_sum;                                            % sampled step limit noise statistics
stats.RT         = RT;                                                      % basis RT data statistics
stats.run_rt     = run_rts;                                                 % run RT stats
stats.dist_rt    = dist_rts;                                                % distamce RT stats
stats.attempt_rt = attempt_rts;                                             % attempt RT stats

end