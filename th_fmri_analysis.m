function th_fmri_analysis(fmri_inps)

% This function evaluates the fMRI data obtained from subjects completing 
% the treasure hunt task in the scanner. 
%
% Inputs
%       fmri_inps:     structure with fields
%        .data_dir:    group-level data directory
%        .func_dir:    functional data subdirectory: participant-level 
%        .glm_dir_l1:  name of the first-level GLM directory
%        .glm_dir_l2:  name of the second-level GLM directory
%        .subs:        BIDS style participant IDs
%        .runs:        participant run numbers
%        .t_in_run:    number of tasks in run 
%        .max_attempt: maximum number of attempts on task  
%        .pomdp:       task pomdp formulation
%        .lumin:       structure with fields - input structure to get stimuli luminance values
%          .stim_dir:  full path: visual stimuli subdirectory
%          .grid_ims:  image names: grid cell position
%          .treas_im:  image name: treasure
%          .obs_im:    image names: observation bars
%          .arrow_ims: image names: decision prompt arrows
%
% Outputs
%       none
% 
% Copyright (C) Lilla Horvath, Dirk Ostwald
% -------------------------------------------------------------------------

% initialization
% -------------------------------------------------------------------------

% invoke SPM12
addpath(genpath([pwd '\SPM12']))

% unpack input structure
data_dir    = fmri_inps.data_dir;
func_dir    = fmri_inps.func_dir;
glm_dir_l1  = fmri_inps.glm_dir_l1;
glm_dir_l2  = fmri_inps.glm_dir_l2;
subs        = fmri_inps.subs;
runs        = fmri_inps.runs;
t_in_run    = fmri_inps.t_in_run; 
max_attempt = fmri_inps.max_attempt;
pomdp       = fmri_inps.pomdp; 
lumin       = fmri_inps.lumin;

% subject-level GLM fMRI analysis
% -------------------------------------------------------------------------

% cycle over subjects
for sj = 1:numel(subs)
    
    % specify input structure
    fl_inps             = [];                                               % initialize
    fl_inps.data_dir    = data_dir;                                         % group-level data directory
    fl_inps.sub         = subs{sj};                                         % subject ID
    fl_inps.func_dir    = func_dir;                                         % subject functional data subdirectory 
    fl_inps.run         = runs{sj};                                         % subject specific run numbers
    fl_inps.t_in_run    = t_in_run;                                         % number of tasks in run
    fl_inps.max_attempt = max_attempt;                                      % maximum number of attempts on task
    fl_inps.pomdp       = pomdp;                                            % task pomdp formulation
    fl_inps.lumin       = lumin;                                            % structure with fields - input structure to get stimuli luminosity values
    fl_inps.glm_dir_l1  = glm_dir_l1;                                       % name of the first-level GLM directory
    
    % formulate, estimate and evaluate the first-level GLM
    spm_fmri_th_glm_first_level(fl_inps)

end

% group-level GLM fMRI analysis
% -------------------------------------------------------------------------

% specify input structure
sl_inps            = [];                                                    % initialize
sl_inps.data_dir   = data_dir;                                              % group-level data directory
sl_inps.subs       = subs;                                                  % subject IDs
sl_inps.glm_dir_l1 = glm_dir_l1;                                            % name of the first-level GLM directory
sl_inps.glm_dir_l2 = glm_dir_l2;                                            % name of the second-level GLM directory

% formulate, estimate and evaluate the second-level GLM
spm_fmri_th_glm_second_level(sl_inps)
    
end

function spm_fmri_th_glm_first_level(fl_inps)

% This function implements the GLM-based fMRI data analysis of a single 
% subject. It performs the following steps:
%
% (1) first-level GLM formulation and GLM visualization as shown in
%     Supplementary Fig. 15 of Horvath et al. (201X) XXX.
% (2) first-level GLM estimation
% (3) first-level GLM evaluation (contrasts of interest)
%
% Inputs
%       fl_inps:       input structure with fields
%        .data_dir:    group-level data root directory
%        .sub:         subject ID
%        .func_dir:    subject functional data subdirectory
%        .run:         subject specific run numbers
%        .t_in_run:    number of tasks in run
%        .max_atempt:  maximum number of attempts on task
%        .pomdp:       task pomdp formulation
%        .lumin:       structure with fields - input structure to get stimuli luminosity values
%          .stim_dir:  visual stimuli subdirectory
%          .grid_ims:  image names: grid cell position
%          .treas_im:  image name: treasure
%          .obs_im:    image names: observation bars
%          .arrow_ims: image names: decision prompt arrows
%        .glm_l1_dir:  name of the first-level GLM directory
%
% Outputs
%       none, writes first-level GLM results to disk      
% 
% Copyright (C) Lilla Horvath, Dirk Ostwald
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%                              initialization 
% -------------------------------------------------------------------------

% unpack input structure
data_dir    = fl_inps.data_dir;
sub         = fl_inps.sub;
func_dir    = fl_inps.func_dir;
run         = fl_inps.run;
t_in_run    = fl_inps.t_in_run;
max_attempt = fl_inps.max_attempt;
pomdp       = fl_inps.pomdp;
lumin       = fl_inps.lumin;
glm_dir_l1  = fl_inps.glm_dir_l1; 

% *************** set up SPM and data management structures ***************

% set SPM defaults
spm('defaults','fmri')
spm_jobman('initcfg');

% get SPM graphics window 
fg = spm_figure('GetWin','Graphics'); 

% inform user
fprintf(['Treasure Hunt ' sub ' GLM\n'])

% full path: subject functional data directory
sub_func_dir = [data_dir, filesep, sub, filesep, func_dir, filesep]; 

% full path: first level GLM target directory that will contain the created .mat and .nii files
tgt_dir_glm_l1 = [data_dir, filesep, sub, filesep glm_dir_l1 filesep];

% delete target directory and its content if pre-existent and create a new one  
if exist(tgt_dir_glm_l1, 'dir')
    rmdir(tgt_dir_glm_l1,'s') 
    mkdir(tgt_dir_glm_l1)
else
    % create target directory if non-existent
    mkdir(tgt_dir_glm_l1)
end

% ******************* convert BIDS data to wide format ********************

% specify input structure 
data_yx         = [];
data_yx.sub_dir = sub_func_dir;                                             % subject data directory
data_yx.sub     = sub;                                                      % subject ID
data_yx.sub_run = run;                                                      % subject run
data_yx.n_task  = t_in_run;                                                 % numbert of tasks in run
data_yx.a_max   = max_attempt;                                              % maximum number of attempts

% create cell arrays y and x containing the run-task-attempt observable and 
% not observable data 
data_yx = th_beh_data_yx(data_yx);

% ************ evaluate trial-by-trial parameters of interest *************

% specify input structure 
params          = [];
params.n_run    = numel(run);                                               % number of runs
params.n_task   = t_in_run;                                                 % number of tasks per run
params.pomdp    = pomdp;                                                    % task pomdp formulation
params.aoi      = @th_agentobj_BB_R_25;                                     % agent of interest: the behaviorally most plausible agent model
params.par_beta = 1.7236;                                                   % group median of the ML softmax parameter (beta) estimate 
params.y        = data_yx.y;                                                % participant observable dataset 

% get the KL divergence and chosen action valence values 
params = th_params_eval(params);

% unpack
KL = params.KL;                                                             % trial-by-trial KL divergence values
Q  = params.Q;                                                              % trial-by-trial chosen action valence values
    
% -------------------------------------------------------------------------
%                            (1) GLM formulation 	
% -------------------------------------------------------------------------

% specify target directory in job structure
job{1}.spm.stats.fmri_spec.dir  = {tgt_dir_glm_l1};

% cycle over runs
for r = 1:max(run)
    
    % ******************* get run specific EPI volumes ********************
    
    % get number of frames 
    N      = [];
    fr_num = [];
    N      = nifti([sub_func_dir, [sub, '_task-th', '_run-0' num2str(r) '_bold.nii']]);
    fr_num = N.dat.dim(4);
     
    % create SPM style filenames 
    f = [];
    f = [sub_func_dir ['snr', sub, '_task-th', '_run-0' num2str(r) '_bold.nii']];
    
    % get list of preprocessed EPI volumes 
    fs = cell(fr_num, 1);
    for fr = 1:fr_num
       fs(fr) = cellstr([f ',' num2str(fr)]);
    end

    % assign list of EPIs to job structure
    job{1}.spm.stats.fmri_spec.sess(r).scans = fs;
    
    % ******************** get run specific regressors ********************
    
    % load the events data file
    doi = [];
    doi = struct2cell(tdfread([sub_func_dir [sub, '_task-th', '_run-0' num2str(r) '_events.tsv']]));
    
    % convert cell to matrix (convert strings to numbers)
    doi_m = [];
    doi_m = zeros(length(doi{1}),numel(doi));
    
    % cycle over the events descriptor fields (row-wise entries of the cell array)
    for k = 1:length(doi)
        % cycle over events
        for l = 1:length(doi{k})
            if ischar(doi{k}(l,:))
                if regexp(doi{k}(l,:), 'n/a')
                    doi_m(l,k) = NaN;                                       % switch n/a strings to numeric NaN variables
                else
                    doi_m(l,k) = str2num(doi{k}(l,:));                      % convert strings to numbers
                end
            else 
                doi_m(l,k) = doi{k}(l);
            end
        end
    end
    
    % unpack
    ons       = doi_m(:,1);                                                 % event onset
    dur       = doi_m(:,2);                                                 % event duration
    ev_type   = doi_m(:,3);                                                 % event type
    trial_num = doi_m(:,4);                                                 % trial number
    pos       = doi_m(:,5);                                                 % agent position
    obs_n     = doi_m(:,6);                                                 % observation bar - north
    obs_e     = doi_m(:,7);                                                 % observation bar - east
    obs_s     = doi_m(:,8);                                                 % observation bar - south
    obs_w     = doi_m(:,9);                                                 % observation bar - west
    act       = doi_m(:,10);                                                % agent action
    tf        = doi_m(:,11);                                                % treasure discovery within attempt
    att_num   = doi_m(:,12);                                                % attempt number
    tk        = doi_m(:,13);                                                % treasure location known from previous attempt
    sli       = doi_m(:,14);                                                % attemtp step limit
    task_num  = doi_m(:,15);                                                % task number
    tgt_1     = doi_m(:,16);                                                % first treasure location
    tgt_2     = doi_m(:,17);                                                % second treasure location
    ops       = doi_m(:,18);                                                % optimal step size 
    
    % concatenate 
    run_events = [];
    run_events = [ons, dur, ev_type, trial_num, pos, obs_n, obs_e, obs_s, obs_w, act, tf, att_num, tk, sli, task_num, tgt_1, tgt_2, ops];
    
    % store run 1 events data for visualization
    if r == 1
        run_events_r1 = run_events;
    end
    
    % basis regressor: full trials
    % ---------------------------------------------------------------------
    
    pos_i = [];                                                             % index marking grid cell position presentation events
    doi   = [];                                                             % all position events
    ioi   = [];                                                             % index of positions of full trials 
    pos_i = find(run_events(:,3)==1);                                 
    doi   = run_events(pos_i,:);
    ioi   = find(~isnan(doi(:,10)));
    
    ftri_ons = [];                                                          % full trial onsets                                            
    ftri_dur = [];                                                          % full trial durations
    ftri_ons = doi(ioi,1);  
    
    % cycle over full trials
    for i = 1:numel(ioi) 
        ftri_dur = [ftri_dur, run_events(pos_i(ioi(i))+3,1)-ftri_ons(i)];   % get durations
    end
    
    % parametric modulators of interest: KL and Q 
    % ---------------------------------------------------------------------
    
    % ------------ partition arrays to run-task-attempt format ------------
    
    % initialize run arrays
    KL_run = [];
    Q_run  = [];
    
    % cycle over tasks
    for t = 1:t_in_run 
        
        % initialize task arrays
        KL_task = [];
        Q_task  = [];
        
        % cycle over attempts
        for a = 1:numel(KL{r}{t})
            % concatenate task values
            KL_task = [KL_task; KL{r}{t}{a}];
            Q_task  = [Q_task; Q{r}{t}{a}];
        end
        
        % concatenate run values
        KL_run = [KL_run; KL_task];
        Q_run  = [Q_run; Q_task];
    end

    % ------------------ get values for full trials only ------------------
    
    % chosen action valence
    Q_run = Q_run(~isnan(Q_run));                                           % get numeric values

    % KL divergence
    tn_run   = [];                                                          % task number to KL values (necessary for task-wise normalization)
    at_found = [];                                                          % attempt treasure found flag (necessary for Q rescaling)
    tn_run   = run_events(pos_i,15); 
    at_found = run_events(pos_i,11);                                     
    
    for i = 1:numel(pos_i)                                                  % cycle over position presentation events
        if isnan(run_events(pos_i(i),10))                                   % if no decision on trial (missed response or last trial on attempt): replace KL value with NaN
            KL_run(i)   = NaN;
            tn_run(i)   = NaN;
            at_found(i) = NaN;
        end
    end
    
    KL_run   = KL_run(~isnan(KL_run));                                      % get numeric values
    tn_run   = tn_run(~isnan(tn_run));                                      % get numeric values
    at_found = at_found(~isnan(at_found));                                  % get numeric values
    
    % multiply chosen action valence by 2 if there is only one treasure to
    % be found in attempt
    Q_run(find(at_found==0)) = Q_run(find(at_found==0))*2;

    %  --------------------- task-wise normalization ----------------------
    
    % initialize task normalized value arrays
    KL_run_norm = [];                                                      
    Q_run_norm  = [];
    
    % cycle over tasks
    for t = 1:4
        
        % get normalized KL 
        nmrtr_KL    = [];
        dnmntr_KL   = [];
        nmrtr_KL    = KL_run(find(tn_run==t))-min(KL_run(find(tn_run==t)));
        dnmntr_KL   = max(KL_run(find(tn_run==t)))-min(KL_run(find(tn_run==t)));
        KL_run_norm = [KL_run_norm; nmrtr_KL/dnmntr_KL];
        
        % get normalized Q 
        nmrtr_Q    = [];
        dnmntr_Q   = [];
        nmrtr_Q    = Q_run(find(tn_run==t))-min(Q_run(find(tn_run==t)));
        dnmntr_Q   = max(Q_run(find(tn_run==t)))-min(Q_run(find(tn_run==t)));
        Q_run_norm = [Q_run_norm; nmrtr_Q/dnmntr_Q];
       
    end
    
    % parametric modulator of no interest: average trial luminance 
    % ---------------------------------------------------------------------
    
    % get stimuli luminance 
    lumin = th_luminance(lumin);
    
    % ------------ compute average luminance of the full trials -----------

    % initialize
    tr_lum_all = [];                                                       
    tr_lum_all = NaN(numel(ioi),1);                                         
    
    % cycle over full trials
    for i = 1:numel(ioi)
        
        % get full trial mean luminance
        tr_lum = zeros(3,1);
        
        % grid position stimuli
        if run_events(pos_i(ioi(i)),5) == 1
            tr_lum(1) = lumin.grid_lum(1);
        elseif run_events(pos_i(ioi(i)),5) == 5
            tr_lum(1) = lumin.grid_lum(2);
        elseif run_events(pos_i(ioi(i)),5) == 21
            tr_lum(1) = lumin.grid_lum(3);
        elseif run_events(pos_i(ioi(i)),5) == 25
            tr_lum(1) = lumin.grid_lum(4);
        elseif ismember(run_events(pos_i(ioi(i)),5), (2:4))
            tr_lum(1) = lumin.grid_lum(5);
        elseif ismember(run_events(pos_i(ioi(i)),5), (6:pomdp.d:16))
            tr_lum(1) = lumin.grid_lum(6);
        elseif ismember(run_events(pos_i(ioi(i)),5), (10:pomdp.d:20))
            tr_lum(1) = lumin.grid_lum(7);
        elseif ismember(run_events(pos_i(ioi(i)),5), (22:24))
            tr_lum(1) = lumin.grid_lum(8);
        else
            tr_lum(1) = lumin.grid_lum(9);
        end
        
        % treasure stimulus
        if ismember(run_events(pos_i(ioi(i)),5), run_events(pos_i(ioi(i)),16:17))
            tr_lum(1) = tr_lum(1)+lumin.treas_lum;
        end
        
        % observation bar stimuli
        % north
        if run_events(pos_i(ioi(i)),6) == 1 
            tr_lum(2) = tr_lum(2) + lumin.obs_lum(2);
        elseif run_events(pos_i(ioi(i)),6) == 0
            tr_lum(2) = tr_lum(2) + lumin.obs_lum(1);
        end
        % east
        if run_events(pos_i(ioi(i)),7) == 1 
            tr_lum(2) = tr_lum(2) + lumin.obs_lum(4);
        elseif run_events(pos_i(ioi(i)),7) == 0
            tr_lum(2) = tr_lum(2) + lumin.obs_lum(3);
        end
        % south
        if run_events(pos_i(ioi(i)),8) == 1 
            tr_lum(2) = tr_lum(2) + lumin.obs_lum(2);
        elseif run_events(pos_i(ioi(i)),8) == 0
            tr_lum(2) = tr_lum(2) + lumin.obs_lum(1);
        end
        % west
        if run_events(pos_i(ioi(i)),9) == 1 
            tr_lum(2) = tr_lum(2) + lumin.obs_lum(4);
        elseif run_events(pos_i(ioi(i)),9) == 0
            tr_lum(2) = tr_lum(2) + lumin.obs_lum(3);
        end
        
        % decision prompt arrow stimuli
        av_act = find(~isnan(run_events(pos_i(ioi(i)),6:9)));
        
        for j = 1:numel(av_act)
            if av_act(j) == 1
                tr_lum(3) = tr_lum(3) + lumin.arrow_lum(1);
            elseif av_act(j) == 2
                tr_lum(3) = tr_lum(3) + lumin.arrow_lum(4);
            elseif av_act(j) == 3
                tr_lum(3) = tr_lum(3) + lumin.arrow_lum(2);
            elseif av_act(j) == 4
                 tr_lum(3) = tr_lum(3) + lumin.arrow_lum(3);
            end
        end
        
        % average trial luminance
        tr_lum_all(i) = mean(tr_lum);
        
    end
    
    % ---------------------- task-wise normalization ----------------------
    
    tr_lum_norm = (tr_lum_all-min(tr_lum_all))./(max(tr_lum_all)-min(tr_lum_all));
    
    % check regressor size equivality 
    % ---------------------------------------------------------------------
    
    % full trial onsets and duration
    if numel(ftri_ons) ~= numel(ftri_dur)
        error('Full trial onset and duration mismatch')
    end
    
    % full trial - normalized Bayesian surprise
    if length(KL_run_norm) ~= numel(ftri_ons)
        error('KL param - Regressor size mismatch')
    end
    
    % full trial - normalized chosen action valence
    if length(Q_run_norm) ~= numel(ftri_ons)
        error('Q param - Regressor size mismatch')
    end
    
    % full trial - normalized luminance
    if length(tr_lum_norm) ~= numel(ftri_ons)
        error('Luminance param - Regressor size mismatch')
    end
    
    % additional regressors 
    % ---------------------------------------------------------------------
    
    % ----------------------- treasure found onsets -----------------------
    
    treas_ons = [];
    for j = 1:numel(pos_i)
        if ismember(run_events(pos_i(j),5), run_events(pos_i(j),16:17))
            treas_ons = [treas_ons, run_events(pos_i(j),1)];
        end
    end
    
    % ------------------------ information onsets -------------------------
    % (post-attempt and post-task)
    
    inf_ons = run_events(find(ismember(run_events(:,3),[6,7])),1);
    
    % ----------------------- fixation cross onsets -----------------------
    % (post-attempt information and post-task information)

    fc_ons = run_events(find(ismember(run_events(:,3),[8,9])),1); 
    
    % assign regressors to job structure
    % ---------------------------------------------------------------------
    
    % --------------- full trials and parametric modulators ---------------
    
    % full trials - boxcar function (onset + duration) 
    job{1}.spm.stats.fmri_spec.sess(r).cond(1).name          = 'ft';        % condition label
    job{1}.spm.stats.fmri_spec.sess(r).cond(1).onset         = ftri_ons;    % condition onsets
    job{1}.spm.stats.fmri_spec.sess(r).cond(1).duration      = ftri_dur;    % condition duration
    job{1}.spm.stats.fmri_spec.sess(r).cond(1).tmod          = 0;           % no temporal modulation    
    
    % parametric modulation - normalized KL divergence
    job{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(1).name  = 'nbs';       % parametric modulation name
    job{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(1).param = KL_run_norm; % parametric modulation values
    job{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(1).poly  = 1;           % polynomial expansion
    job{1}.spm.stats.fmri_spec.sess(r).cond(1).orth          = 0;           % no orthogonalization of parametric modulation regressors 
    
    % parametric modulation - normalized chosen action Q
    job{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(2).name  = 'nq';        % parametric modulation name
    job{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(2).param = Q_run_norm;  % parametric modulation values
    job{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(2).poly  = 1;           % polynomial expansion
    job{1}.spm.stats.fmri_spec.sess(r).cond(1).orth          = 0;           % no orthogonalization of parametric modulation regressors 
    
    % parametric modulation - average luminance
    job{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(3).name  = 'lum';       % parametric modulation name
    job{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(3).param = tr_lum_norm; % parametric modulation values
    job{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(3).poly  = 1;           % polynomial expansion
    job{1}.spm.stats.fmri_spec.sess(r).cond(1).orth          = 0;           % no orthogonalization of parametric modulation regressors 
    
    % ------------------------ treasure discovery -------------------------
    
    % treasure found at grid cell position - stick function (onset) 
    job{1}.spm.stats.fmri_spec.sess(r).cond(2).name     = 'treas';          % condition label
    job{1}.spm.stats.fmri_spec.sess(r).cond(2).onset    = treas_ons;        % condition onsets
    job{1}.spm.stats.fmri_spec.sess(r).cond(2).duration = 0;                % condition duration
    job{1}.spm.stats.fmri_spec.sess(r).cond(2).tmod     = 0;                % no temporal modulation  
    job{1}.spm.stats.fmri_spec.sess(r).cond(2).orth     = 0;    

    % no parametric modulation
    job{1}.spm.stats.fmri_spec.sess(r).cond(2).pmod     = struct('name', {}, 'param', {}, 'poly', {}); 
    
    % --------------------- information presentation ----------------------
    
    % post-attempt and post-task information - stick function (onset) 
    job{1}.spm.stats.fmri_spec.sess(r).cond(3).name     = 'inf';            % condition label
    job{1}.spm.stats.fmri_spec.sess(r).cond(3).onset    = inf_ons;          % condition onsets
    job{1}.spm.stats.fmri_spec.sess(r).cond(3).duration = 0;                % condition duration
    job{1}.spm.stats.fmri_spec.sess(r).cond(3).tmod     = 0;                % no temporal modulation  
    job{1}.spm.stats.fmri_spec.sess(r).cond(3).orth     = 0;
    
    % no parametric modulation
    job{1}.spm.stats.fmri_spec.sess(r).cond(3).pmod     = struct('name', {}, 'param', {}, 'poly', {});  
    
    % --------------------------- fixation cross --------------------------
    
    % fixation cross after post-attempt and post-task information - stick function (onset) 
    job{1}.spm.stats.fmri_spec.sess(r).cond(4).name     = 'fc';             % condition label
    job{1}.spm.stats.fmri_spec.sess(r).cond(4).onset    = fc_ons;           % condition onsets
    job{1}.spm.stats.fmri_spec.sess(r).cond(4).duration = 0;                % condition duration
    job{1}.spm.stats.fmri_spec.sess(r).cond(4).tmod     = 0;                % no temporal modulation
    job{1}.spm.stats.fmri_spec.sess(r).cond(4).orth     = 0;
    
    % no parametric modulation
    job{1}.spm.stats.fmri_spec.sess(r).cond(4).pmod     = struct('name', {}, 'param', {}, 'poly', {});  
       
    % ----------------------- nuissance regressors ------------------------
    
    % movement parameters
    job{1}.spm.stats.fmri_spec.sess(r).hpf       = 128;                     % high-pass filter cut-off for cosine-set based GLM data filtering = % 1/128 Hz = 0.0078 Hz (default)
    job{1}.spm.stats.fmri_spec.sess(r).multi     = {''};                    % no multiple condition specification
    job{1}.spm.stats.fmri_spec.sess(r).regress   = struct('name', {}, 'val', {}); % no additional regressors
    job{1}.spm.stats.fmri_spec.sess(r).multi_reg = {[sub_func_dir, ['rp_', sub, '_task-th', '_run-0' num2str(r) '_bold.txt']]}; % get realignment parameters
 
end
        
% ---------- additional run and condition independent parameters ----------

job{1}.spm.stats.fmri_spec.timing.units     = 'secs';                       % temporal units
job{1}.spm.stats.fmri_spec.timing.RT        = 2;                            % repetition time (TR)
job{1}.spm.stats.fmri_spec.timing.fmri_t    = 16;                           % microtime resolution for HRF convolution (default)
job{1}.spm.stats.fmri_spec.timing.fmri_t0   = 1;                            % microtime onset for HRF convolution (default)
job{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];                        % no modelling of HRF derivatives -> temporal and dispertion derivatives could go here
job{1}.spm.stats.fmri_spec.volt             = 1;                            % no interaction modeling for factorial designs 
job{1}.spm.stats.fmri_spec.global           = 'None';                       % no global normalization
job{1}.spm.stats.fmri_spec.mask             = {''};                         % no mask
job{1}.spm.stats.fmri_spec.cvi              = 'AR(1)';                      % AR(1) model-based covariance components for pre-whitening
job{1}.spm.stats.fmri_spec.fact             = struct('name', {}, 'levels', {}); % no factorial design

% ********************** formulate GLM design matrix **********************

% inform user
fprintf(['Treasure Hunt ' sub ' Formalizing GLM\n'])

% run job
spm_jobman('run',job);

% clear job 
clear job

% ***** plot run 1 regressors in time-series format (Suppl. Fig. 15a) *****

% initialization
% -------------------------------------------------------------------------

% -------------------------- get GLM design data --------------------------

% load SPM.mat file
load(fullfile(tgt_dir_glm_l1,'SPM.mat'))

% unpack
X               = SPM.xX.X;                                                 % complete design matrix (HRF convolved regressor values)                                              
nscan           = [1 cumsum(SPM.nscan)];                                    % run specific scan indices
U_r1            = SPM.Sess(1).U;                                            % run 1 regressor structure specification
dt              = U_r1.dt;                                                  % time bin length
rep_t           = SPM.xY.RT;                                                % TR 
unc_reg_r1_1    = U_r1(1).u;                                                % unconvolved regressor: full trial, KL, Q, luminance (boxcar)
unc_reg_r1_2    = full(U_r1(2).u);                                          % unconvolved regressor: treasure (stick)
unc_reg_r1_3    = full(U_r1(3).u);                                          % unconvolved regressor: information (stick)
unc_reg_r1_4    = full(U_r1(4).u);                                          % unconvolved regressor: fixation cross (stick)
fMRI_T          = SPM.xBF.T;                                                % number of subdivisions of TR 
fMRI_T0         = SPM.xBF.T0;                                               % first time bin

% regressor downsampling - full trial, KL, Q, luminance
d_unc_reg_r1_1 = full(unc_reg_r1_1((0:(nscan(2) - 1))*fMRI_T + fMRI_T0 + 32,:));

% microtime vector
t_hr = -32*dt:dt:dt*(length(unc_reg_r1_1)-33);

% regressor downsampling - treasure, information, fixation cross
d_unc_reg_r1_2 = zeros(nscan(2),1);
d_unc_reg_r1_3 = zeros(nscan(2),1);
d_unc_reg_r1_4 = zeros(nscan(2),1);
d_unc_reg_r1_2(round(t_hr(find(unc_reg_r1_2))/rep_t)) = 0.17;
d_unc_reg_r1_3(round(t_hr(find(unc_reg_r1_3))/rep_t)) = 0.17;
d_unc_reg_r1_4(round(t_hr(find(unc_reg_r1_4))/rep_t)) = 0.17;

% concatenate downsampled unconvolved regressors (run 1)
d_unc_reg_r1 = [d_unc_reg_r1_1 d_unc_reg_r1_2 d_unc_reg_r1_3 d_unc_reg_r1_4]; 

% regressor labels
reglab = {'Full trial', 'Normalized Bayesian surprise', 'Normalized chosen valence', 'Normalized luminance', 'Treasure found', 'Information display', 'Fixation cross'} ; % regressor labels
nreg   = numel(reglab);                                                    % number of regressors

% ------------------ map task and attempt starts to TRs -------------------

% identify task start TRs
task_vols = [ round(run_events_r1(find(run_events_r1(:,15)==1,1),1)/2),...
              round(run_events_r1(find(run_events_r1(:,15)==2,1),1)/2),...
              round(run_events_r1(find(run_events_r1(:,15)==3,1),1)/2),...
              round(run_events_r1(find(run_events_r1(:,15)==4,1),1)/2) ];            

% identify attempt start TRs
att_vols = [];
for t = 1:t_in_run                                                          % cycle over tasks
    doi = [];
    doi = run_events_r1(find(run_events_r1(:,15)==t),:); 
    for a = 1:numel(KL{1}{t})                                               % cycle over attempts
        att_vols = [att_vols, round(doi(find(doi(:,12)==a,1),1)/2)];
    end
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% TO DO: KEEP ONE OF THE FOLLOWING TWO VIS. VERSIONS!!! 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ***
% (1) LARGE: Dirk's print version
% ***

% --------------------------- initialize figure ---------------------------

h = figure;
set(h, 'Color', [1, 1, 1]);
cols = ([255, 44, 44; 153, 153, 153])./255;                                 % plot colors

% visualization
% ------------------------------------------------------------------------

% cycle over regressors
for reg = 1:nreg
    
    % ---------------------------- get subplot ----------------------------
    
    subplot(nreg,1,reg)
    
    % convolved regressor 
    plot(nscan(1):nscan(2), X(nscan(1):nscan(2),reg), 'LineWidth', 2, 'Color', 'k') 
    hold on
    
    % unconvolved regressor
    plot(nscan(1):nscan(2), d_unc_reg_r1(:,reg), 'LineWidth', 1, 'Color', cols(2,:)) 
    hold on 
    
    % task and attempt start markers 
    plot(att_vols, ones(length(att_vols))*min(X(nscan(1):nscan(2),reg)), 's', 'MarkerFaceColor', cols(1,:), 'MarkerEdgeColor', cols(1,:), 'MarkerSize', 8)
    hold on
    plot(task_vols, ones(length(task_vols))*min(X(nscan(1):nscan(2),reg)), 'o', 'MarkerFaceColor', cols(1,:), 'MarkerEdgeColor', cols(1,:), 'MarkerSize', 10)
    
    % additional subplot formatting
    xlim([nscan(1) nscan(2)])
    if reg > size(unc_reg_r1_1,2)
        set(gca, 'ytick', 0:0.15:0.3, 'FontName', 'Arial', 'FontSize', 10)
         ylim([-inf 0.34])
    else
        set(gca, 'ytick', 0:0.5:1, 'FontName', 'Arial', 'FontSize', 10)
        ylim([-inf 1.2])
    end
    title(reglab{reg}, 'fontweight','normal', 'FontSize', 12) 
    box off
    if reg == nreg
        xlabel('TR', 'FontName', 'Arial', 'FontSize', 12)  
    end
    
end

% save figure 
set(h, 'Position', get(0,'Screensize'), 'PaperPositionMode', 'auto');
print(h, '-dbmp', '-r300', fullfile(tgt_dir_glm_l1, 'SFig15a_v1'))
close(h)

% ***
% (2) SMALL: Lilla's export_fig version
% ***

% --------------------------- initialize figure ---------------------------

% invoke Yair Altman's export_fig toolbox (usage: save figures with minimal white space)
addpath(genpath([pwd '\export_fig']))

% set figure position and background color
h = figure;
set(h, 'Color', [1, 1, 1], 'Units', 'centimeters');
fig_pos     = get(h,'position');
fig_pos_adj = [fig_pos(1:2), 13, 13/fig_pos(3)*fig_pos(4)];
set(h, 'position', fig_pos_adj)

% plot colors
cols = ([255, 44, 44; 153, 153, 153])./255;                                 

% visualization
% ------------------------------------------------------------------------

% cycle over regressors
for reg = 1:nreg
    
    % ---------------------------- get subplot ----------------------------
    
    subplot(nreg,1,reg)
    
    % convolved regressor 
    plot(nscan(1):nscan(2), X(nscan(1):nscan(2),reg), 'LineWidth', 0.5, 'Color', 'k') 
    hold on
    
    % unconvolved regressor
    plot(nscan(1):nscan(2), d_unc_reg_r1(:,reg), 'LineWidth', 0.5, 'Color', cols(2,:)) 
    hold on 
    
    % task and attempt start markers 
    plot(att_vols, ones(length(att_vols))*min(X(nscan(1):nscan(2),reg)), 's', 'MarkerFaceColor', cols(1,:), 'MarkerEdgeColor', cols(1,:), 'MarkerSize', 3)
    hold on
    plot(task_vols, ones(length(task_vols))*min(X(nscan(1):nscan(2),reg)), 'o', 'MarkerFaceColor', cols(1,:), 'MarkerEdgeColor', cols(1,:), 'MarkerSize', 4)
    
    % additional subplot formatting
    xlim([nscan(1) nscan(2)])
    if reg > size(unc_reg_r1_1,2)
        ylim([-inf 0.34])
        set(gca, 'ytick', 0:0.15:0.3)
    else
        ylim([-inf 1.2])
        set(gca, 'ytick', 0:0.5:1)
    end
    set(gca, 'FontName', 'Arial', 'FontSize', 5) 
    title(reglab{reg}, 'fontweight','normal', 'FontSize', 6) 
    box off
    if reg == nreg
        xlabel('TR', 'FontName', 'Arial', 'FontSize', 6)  
    end
    
end

% save figure 
sfig15a_name = fullfile(tgt_dir_glm_l1, ['SFig15a_v2', '.pdf']);
export_fig(h, sfig15a_name);
close(h)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% END OF TO DO
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% *** plot GLM design matrix for the complete dataset (Suppl. Fig. 15b) ***

% save SPM window - Supplementary Fig. 15b
saveas(fg, fullfile(tgt_dir_glm_l1, ['SFig15b', '.pdf']));

% -------------------------------------------------------------------------
%                             (2) GLM estimation 
% -------------------------------------------------------------------------

% *********************** job structure allocation ************************

% GLM to estimate
job{1}.spm.stats.fmri_est.spmmat = {fullfile(tgt_dir_glm_l1,'SPM.mat')}; 

% classical ReML estimation
job{1}.spm.stats.fmri_est.method.Classical = 1; 

% ***************************** estimate GLM ******************************

% inform user
fprintf(['Treasure Hunt ' sub ' Estimating GLM\n'])

% run job
spm_jobman('run', job);

% clear job 
clear job

% -------------------------------------------------------------------------
%                            (3) GLM evaluation 
% -------------------------------------------------------------------------

% *********************** job structure allocation ************************

% GLM to evaluate
job{1}.spm.stats{1}.con.spmmat = {fullfile(tgt_dir_glm_l1,'SPM.mat')}; 

% contrast of interest names
connames = { 'FT',          ...                                             % (1) pos main effect: full trial 
             'NBS',         ...                                             % (2) pos main effect: normalized Bayesian surprise
             'NQ',          ...                                             % (3) pos main effect: normalized chosen action valence
             'LUM',         ...                                             % (4) pos main effect: normalized luminance
             'TREAS',       ...                                             % (5) pos main effect: treasure found
             'INF',         ...                                             % (6) pos main effect: post-attempt and post-task information
             'FC',          ...                                             % (7) pos main effect: fixation cross block
             'nFT',         ...                                             % (8) neg main effect: normalized Bayesian surprise
             'nNBS',        ...                                             % (9) neg main effect: normalized Bayesian surprise
             'nNQ',         ...                                             % (10) neg main effect: normalized chosen action valence
             'nLUM',        ...                                             % (11) neg main effect: normalized luminance
             'nTREAS',      ...                                             % (12) neg main effect: treasure found
             'nINF',        ...                                             % (13) neg main effect: post-attempt and post-task information
             'nFC',         ...                                             % (14) neg main effect: fixation cross block
             'FT_l_NBS',    ...                                             % (15) directed contrast: full trial > normalized Bayesian surprise
             'FT_s_NBS',    ...                                             % (16) directed contrast: full trial < normalized Bayesian surprise
             'NBS_l_NQ',    ...                                             % (17) directed contrast: normalized Bayesian surprise > normalized chosen action valence
             'NBS_s_NQ',    ...                                             % (18) directed contrast: normalized Bayesian surprise < normalized chosen action valence
             'FT_l_NQ',     ...                                             % (19) directed contrast: full trial > normalized chosen action valence
             'FT_s_NQ',     ...                                             % (20) directed contrast: full trial < normalized chosen action valence
             'FT_l_NBS_NQ', ...                                             % (21) directed contrast: full trial > normalized Bayesian surprise + normalized chosen action valence
             'FT_s_NBS_NQ', ...                                             % (22) directed contrast: full trial < normalized Bayesian surprise + normalized chosen action valence
             'NBS_l_FT_NQ', ...                                             % (23) directed contrast: normalized Bayesian surprise > full trial + normalized chosen action valence
             'NBS_s_FT_NQ', ...                                             % (24) directed contrast: normalized Bayesian surprise < full trial + normalized chosen action valence
             'NQ_l_FT_NBS', ...                                             % (25) directed contrast: normalized chosen action valence > full trial + normalized Bayesian surprise
             'NQ_s_FT_NBS'      };                                          % (26) directed contrast: normalized chosen action valence < full trial + normalized Bayesian surprise

% contrast vectors
conweights = [num2cell(eye(nreg),1),    ...
              num2cell(eye(nreg)*-1,1), ... 
              { [1 -1 0 0 0 0 0],       ... 
                [-1 1 0 0 0 0 0],       ...
                [0 1 -1 0 0 0 0],       ...
                [0 -1 1 0 0 0 0],       ...
                [1 0 -1 0 0 0 0],       ...
                [-1 0 1 0 0 0 0],       ...
                [2 -1 -1 0 0 0 0],      ...
                [-2 1 1 0 0 0 0],       ...
                [-1 2 -1 0 0 0 0],      ...
                [1 -2 1 0 0 0 0],       ...
                [-1 -1 2 0 0 0 0],      ...
                [1 1 -2 0 0 0 0] }          ];
                               
% cycle over contrast specifications
for c = 1:numel(connames)
 
    % allocate t-contrast structure
    job{1}.spm.stats{1}.con.consess{c}.tcon.name    = connames{c};       
    job{1}.spm.stats{1}.con.consess{c}.tcon.weights = conweights{c};         
    job{1}.spm.stats{1}.con.consess{c}.tcon.sessrep = 'repl';             

    % delete existing contrasts
    job{1}.spm.stats{1}.con.delete = 1;   
     
end 

% ************************ GLM parameter inference ************************

% inform user
fprintf(['Treasure Hunt ' sub ' Evaluating GLM\n'])

% run job
spm_jobman('run', job);

% clear job 
clear job

end

function spm_fmri_th_glm_second_level(sl_inps)

% This function implements the GLM-based fMRI data analysis of the group of 
% subjects: using the summary statistics approach the group-level data is 
% evaluated based on one-sample t-tests of the single subject contrasted 
% beta estimates. It performs the following steps:
%
% (1) second-level GLM formulation
% (2) second-level GLM estimation
% (3) second-level GLM evaluation (contrasts of interest)
% (4) results visualization
%
% Inputs
%       sl_inps:      input structure with fields
%        .data_dir:   group-level data directory
%        .subs:       subject IDs
%        .glm_dir_l1: name of the first-level GLM directory
%        .glm_dir_l2: name of the second-level GLM directory
%
% Outputs
%       none, writes second-level GLM results to disk      
% 
% Copyright (C) Dirk Ostwald, Lilla Horvath
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%                             initialization 
% -------------------------------------------------------------------------

% unpack input structure
data_dir   = sl_inps.data_dir;
subs       = sl_inps.subs;
glm_dir_l1 = sl_inps.glm_dir_l1;
glm_dir_l2 = sl_inps.glm_dir_l2;

% *************** set up SPM and data management structures ***************

% set SPM defaults
spm('defaults','fmri')
spm_jobman('initcfg');

% full path: second level GLM target directory that will contain the created job.mat and SPM.mat files 
tgt_dir = [data_dir, filesep glm_dir_l2 filesep];

% delete target directory and its content if pre-existent and create a new one
if exist(tgt_dir, 'dir')
    rmdir(tgt_dir,'s')
    mkdir(tgt_dir)
else
    % create target directory if non-existent
    mkdir(tgt_dir)
end

% *********************** define GLM model contrasts **********************

% contrast of interest names
connames = { 'FT',          ...                                             % (1) pos main effect: full trial 
             'NBS',         ...                                             % (2) pos main effect: normalized Bayesian surprise
             'NQ',          ...                                             % (3) pos main effect: normalized chosen action valence
             'LUM',         ...                                             % (4) pos main effect: normalized luminance
             'TREAS',       ...                                             % (5) pos main effect: treasure found
             'INF',         ...                                             % (6) pos main effect: post-attempt and post-task information
             'FC',          ...                                             % (7) pos main effect: fixation cross block
             'nFT',         ...                                             % (8) neg main effect: normalized Bayesian surprise
             'nNBS',        ...                                             % (9) neg main effect: normalized Bayesian surprise
             'nNQ',         ...                                             % (10) neg main effect: normalized chosen action valence
             'nLUM',        ...                                             % (11) neg main effect: normalized luminance
             'nTREAS',      ...                                             % (12) neg main effect: treasure found
             'nINF',        ...                                             % (13) neg main effect: post-attempt and post-task information
             'nFC',         ...                                             % (14) neg main effect: fixation cross block
             'FT_l_NBS',    ...                                             % (15) directed contrast: full trial > normalized Bayesian surprise
             'FT_s_NBS',    ...                                             % (16) directed contrast: full trial < normalized Bayesian surprise
             'NBS_l_NQ',    ...                                             % (17) directed contrast: normalized Bayesian surprise > normalized chosen action valence
             'NBS_s_NQ',    ...                                             % (18) directed contrast: normalized Bayesian surprise < normalized chosen action valence
             'FT_l_NQ',     ...                                             % (19) directed contrast: full trial > normalized chosen action valence
             'FT_s_NQ',     ...                                             % (20) directed contrast: full trial < normalized chosen action valence
             'FT_l_NBS_NQ', ...                                             % (21) directed contrast: full trial > normalized Bayesian surprise + normalized chosen action valence
             'FT_s_NBS_NQ', ...                                             % (22) directed contrast: full trial < normalized Bayesian surprise + normalized chosen action valence
             'NBS_l_FT_NQ', ...                                             % (23) directed contrast: normalized Bayesian surprise > full trial + normalized chosen action valence
             'NBS_s_FT_NQ', ...                                             % (24) directed contrast: normalized Bayesian surprise < full trial + normalized chosen action valence
             'NQ_l_FT_NBS', ...                                             % (25) directed contrast: normalized chosen action valence > full trial + normalized Bayesian surprise
             'NQ_s_FT_NBS'      };                                          % (26) directed contrast: normalized chosen action valence < full trial + normalized Bayesian surprise

% cycle over contrasts
% -------------------------------------------------------------------------
for con = 1:numel(connames)
    
    % ---------------------------------------------------------------------
    %                       (1) GLM formulation 
    % ---------------------------------------------------------------------
    
    % ******************* specify data for GLM analysis *******************
    
    % initialize filename array
    fNames = [];
    
    % cycle over subjects
    for sj = 1:numel(subs)
        
        % full path: first level GLM data directory
        sub_glm_dir = fullfile(data_dir, subs{sj}, glm_dir_l1);
        
        % create SPM style filenames
        if con < 10
            f = ['con_000' num2str(con) '.nii'];
        else
            f = ['con_00'  num2str(con) '.nii'];
        end
        
        fs = cellstr([repmat([sub_glm_dir filesep], 1, 1) f repmat(',1', 1, 1)]);
        
        % concatenate over subjects
        fNames = [fNames; fs];
        
    end
    
    % contrast subdirectory for each second level analysis
    con_dir = fullfile(tgt_dir, ['Con' num2str(con) '_' connames{con}]); 
    
    % delete contrast directory and its content if pre-existent and create a new one
    if exist(con_dir, 'dir')
        rmdir(con_dir,'s')
        mkdir(con_dir)
    else
        % create contrast subdirectory if non-existent
        mkdir(con_dir)
    end
    
    % ********************** job structure allocation *********************
    job{1}.spm.stats.factorial_design.dir                     = {con_dir};
    job{1}.spm.stats.factorial_design.des.t1.scans            = cellstr(fNames);
    job{1}.spm.stats.factorial_design.cov                     = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    job{1}.spm.stats.factorial_design.masking.tm.tm_none      = 1;
    job{1}.spm.stats.factorial_design.masking.im              = 1;
    job{1}.spm.stats.factorial_design.masking.em              = {''};
    job{1}.spm.stats.factorial_design.globalc.g_omit          = 1;
    job{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no  = 1;
    job{1}.spm.stats.factorial_design.globalm.glonorm         = 1;
    
    % *********** formulate one sample t-test GLM design matrix ***********
    
    % run job
    spm_jobman('run', job)
    
    % clear job
    clear job
    
    % ---------------------------------------------------------------------
    %                        (2) GLM estimation 
    % ---------------------------------------------------------------------
    
    % ********************** job structure allocation *********************

    % GLM to estimate
    job{1}.stats{1}.fmri_est.spmmat = {fullfile(con_dir, 'SPM.mat')}; 

    % *************************** estimate GLM ****************************
    
    % run job
    spm_jobman('run', job);
    
    % clear job
    clear job
    
    % ---------------------------------------------------------------------
    %                         (3) GLM evaluation 
    % ---------------------------------------------------------------------
   
    % ********************** job structure allocation *********************
    
    % GLM to evaluate
    job{1}.spm.stats{1}.con.spmmat = {fullfile(con_dir, 'SPM.mat')}; 

    % allocate t-contrast structure
    job{1}.spm.stats{1}.con.consess{1}.tcon.name    = connames{con};       
    job{1}.spm.stats{1}.con.consess{1}.tcon.weights = 1;         
       
    % delete existing contrasts
    job{1}.spm.stats{1}.con.delete = 1;   
    
    % ********************** GLM parameter inference **********************

    % run job 
    spm_jobman('run', job);
    
    % clear job
    clear job

    % ---------------------------------------------------------------------
    %                       (4) results visualization 
    % ---------------------------------------------------------------------
    
    % ********************** job structure allocation *********************
    
    job{1}.spm.stats.results.spmmat              = {fullfile(con_dir, 'SPM.mat')}; % GLM to evaluate
    job{1}.spm.stats.results.conspec.titlestr    = connames{con};           % contrast name
    job{1}.spm.stats.results.conspec.contrasts   = 1;                       % contrast index
    job{1}.spm.stats.results.conspec.threshdesc  = 'none';                  % multiple testing control
    job{1}.spm.stats.results.conspec.thresh      = 0.001;                   % p-value threshold
    job{1}.spm.stats.results.conspec.extent      = 0;                       % cluster extent treshold
    job{1}.spm.stats.results.conspec.conjunction = 1;                       % conjunction analysis setting
    job{1}.spm.stats.results.conspec.mask.none   = 1;                       % no mask
    job{1}.spm.stats.results.units               = 1;                       % units        
    job{1}.spm.stats.results.print               = 'fig';                   % figure type
    job{1}.spm.stats.results.write.none          = 1;                       % no result writing

    % ************************* visualize results *************************
    
    % run job
    spm_jobman('run', job);
    
    % clear job
    clear job
    
    close all
end
end

