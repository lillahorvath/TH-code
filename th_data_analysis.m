function th_data_analysis

% This function implements all analysis reported in Horvath et al. (201X) 
% 'The XXX'. These are the following:
%
% (1) treasure hunt experimental task formulation as PoMDP
% (2) behavioral data analysis - participants and agents
% (3) fMRI preprocessing
% (4) fMRI data analysis
% 
% Inputs
%       none
%
% Outputs
%       none
% 
% Copyright (C) Lilla Horvath, Dirk Ostwald
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%                             initialization 
% -------------------------------------------------------------------------

% initialize environment
clc
clear all
close all 
dbstop if error

% ************************* directory management **************************

% main data directory: contains the BIDS formatted data and stored outputs of the 
% analyses
data_dir = 'C:\Users\Lilla\Documents\Treasure_hunt\Data\ds_xxx';

% data subdirectories and stored file names
br_dir     = 'beh_res';                                                     % behavioral results subdirectory: group-level 
func_dir   = 'func';                                                        % functional data subdirectory: participant-level 
pre_dir    = 'pre';                                                         % fmri preprocessing subdirectory: participant-level 
glm_dir_l1 = 'glm_l1';                                                      % fMRI results subdirectory: participant-level 
glm_dir_l2 = 'glm_l2';                                                      % fMRI results subdirectory: group-level 
pomdp_dir  = 'pomdp';                                                       % pomdp task formulation results subdirectory
stim_dir   = 'stimuli';                                                     % presented visual stimuli subdirectory
p_brn      = 'th_beh_res_p.mat';                                            % name of the behavioral results file - participants
a_brn      = 'th_beh_res_a.mat';                                            % name of the behavioral results file - agents

% participant IDs
subs = { 'sub-01', ...
         'sub-02', ...
         'sub-03', ...
         'sub-04', ...
         'sub-05', ...
         'sub-06', ...
         'sub-07', ...
         'sub-08', ...
         'sub-09', ...
         'sub-10', ...
         'sub-11', ...
         'sub-12', ...
         'sub-13', ...
         'sub-14', ...
         'sub-15', ...
         'sub-16', ...
         'sub-17', ...
         'sub-18', ...
         'sub-19'      };

% participant-specific run numbering
runs = { [1 2 3 4], ...
         [1 2 3 4], ...
         [1 2 3 4], ...
         [1 2 3 4], ...
         [1 2 3 4], ...
         [1 2 3 4], ...
         [1 2 3 4], ...
         [1 2 3  ], ...
         [1 2 3 4], ...
         [1 2 3 4], ...
         [1 2 3 4], ...
         [1 2 3 4], ...
         [1 2 3 4], ...
         [1 2 3 4], ...
         [1 2 3 4], ...
         [1 2 3 4], ...
         [1 2 3 4], ...
         [1 2 3  ], ...
         [1 2 3 4]      };
                            
% *************************** agents of interest **************************

% agent model handles
agents = { @th_agentobj_BF_0,       ...                                     % agent BF-0
           @th_agentobj_BF_R,       ...                                     % agent BF-R
           @th_agentobj_BB_R_1,     ...                                     % agent BB-R-1
           @th_agentobj_BB_R_2,     ...                                     % agent BB-R-2
           @th_agentobj_BB_R_25,    ...                                     % agent BB-R-25
           @th_agentobj_BB_R_25_2,  ...                                     % agent BB-R-25-2
           @th_agentobj_BB_I,       ...                                     % agent BB-I
           @th_agentobj_BB_H_C,     ...                                     % agent BB-H-C
           @th_agentobj_BB_H_E          };                                  % agent BB-H-E 
                                     
% ************************** experiment constants *************************

% basic experiment and task parameters
max_run     = 4;                                                            % maximum number of runs
t_in_run    = 4;                                                            % number of tasks in run
max_attempt = 3;                                                            % maximum number of attempts
max_opt_dis = 12;                                                           % maximum length of shortest path (= optimal step size)
sln_set     = 5;                                                            % set size of the step limit noise 
dim         = 5;                                                            % grid-world dimensionality

% visual stimuli filenames
grid_fname  = { 'grid_upper_corner_left.jpg',  ...                          % grid cell position
                'grid_upper_corner_right.jpg', ...
                'grid_lower_corner_left.jpg',  ...
                'grid_lower_corner_right.jpg', ...
                'grid_side_upper.jpg',         ...
                'grid_side_left.jpg',          ...
                'grid_side_right.jpg',         ...
                'grid_side_lower.jpg',         ...
                'grid_middle.jpg'                  };                                 
obs_fname   = { 'ob_dark_hor.jpg',             ...                          % observation bars  
                'ob_light_hor.jpg',            ...
                'ob_dark_ver.jpg',             ...
                'ob_light_ver.jpg'                 };               
arrow_fname = { 'dpa_up.jpg',                  ...                          % decision prompt arrows
                'dpa_down.jpg',                ...
                'dpa_left.jpg',                ...                    
                'dpa_right.jpg'                    };
treas_fname =   'treasure.jpg';                                             % treasure     

% ************************** analysis parameters **************************

% agent simulation parameters
nrep       = 1e0;                                                           % number of simulation repeats
sim_beta   = 1e2;                                                           % softmax beta parameter (all agents except BF-0)
sim_lambda = 0.7;                                                           % linear convex combination lambda parameter (BB-H-C)
sim_os_sl  = [1, 0.4];                                                      % exponential decay offset and slope parameter (BB-H-E)
toi        = [8 25];                                                        % target locations of interest
stoi       = 0;                                                             % flag: evaluate agent only on task with target locations of interest  

% agent maximum log likelihood evaluation - optimization parameters
beta_min    = 1e-5;                                                         % softmax beta lower interval endpoint (all agents except BF-0)     
beta_max    = 1e1;                                                          % softmax beta upper interval endpoint (all agents except BF-0)
beta_init   = 5e0;                                                          % softmax beta initial point (all agents except BF-0)
lambda_min  = 1e-5;                                                         % linear convex combination lambda lower interval endpoint (BB-H-C)   
lambda_max  = 1e0;                                                          % linear convex combination lambda upper interval endpoint (BB-H-C)
lambda_init = 0.5;                                                          % linear convex combination lambda initial point (BB-H-C)
os_min      = 1e-5;                                                         % exponential decay offset lower interval endpoint (BB-H-E)
os_max      = 1e0;                                                          % exponential decay offset upper interval endpoint (BB-H-E)
os_init     = 0.5;                                                          % exponential decay offset initial point (BB-H-E)  
sl_min      = 1e-5;                                                         % exponential decay slope lower interval endpoint (BB-H-E)
sl_max      = 1e1;                                                          % exponential decay slope upper interval endpoint (BB-H-E)
sl_init     = 5;                                                            % exponential decay slope initial point (BB-H-E)

% visualization parameters
plot_pomdp = 1;                                                             % flag: visualize pomdp task formulation

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% TO DO: check th_agent_fit before turning on plot_ll
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_ll    = 0;                                                             % flag: plot log likelihood function 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% END OF TO DO
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% -------------------------------------------------------------------------
%                        (1) PoMDP task formulation 
% -------------------------------------------------------------------------

% specify input structure
pomdp   = [];                                                               % initialize input structure 
pomdp.d = dim;                                                              % grid-world dimensionality   

% formulate task as PoMDP
pomdp = th_pomdp(pomdp, data_dir, pomdp_dir, plot_pomdp);   

% -------------------------------------------------------------------------
%                       (2) behavioral data analysis 
% -------------------------------------------------------------------------

% specify input structure
beh_inps             = [];                                                  % initialize input structure
beh_inps.data_dir    = data_dir;                                            % group-level (=main) data directory
beh_inps.func_dir    = func_dir;                                            % functional data subdirectory: participant-level
beh_inps.br_dir      = br_dir;                                              % behavioral results subdirectory: group-level 
beh_inps.p_brn       = p_brn;                                               % name of the behavioral results file - participants
beh_inps.a_brn       = a_brn;                                               % name of the behavioral results file - agents
beh_inps.subs        = subs;                                                % BIDS style participant IDs
beh_inps.runs        = runs;                                                % participant run numbers
beh_inps.max_run     = max_run;                                             % maximum number of experimental runs 
beh_inps.t_in_run    = t_in_run;                                            % number of tasks in run 
beh_inps.max_attempt = max_attempt;                                         % maximum number of attempts on task
beh_inps.max_opt_dis = max_opt_dis;                                         % maximum length of shortest path (= optimal step size) 
beh_inps.sln_set     = sln_set;                                             % set size of the step limit noise 
beh_inps.dim         = dim;                                                 % grid-world dimensionality         
beh_inps.pomdp       = pomdp;                                               % task pomdp formulation
beh_inps.agents      = agents;                                              % agents of interest                            
beh_inps.nrep        = nrep;                                                % number of simulation repeats
beh_inps.sim_beta    = sim_beta;                                            % softmax beta parameter for agent simulations (all agents except BF-0)
beh_inps.sim_lambda  = sim_lambda;                                          % linear convex combination lambda parameter for agent simulation (BB-H-C)
beh_inps.sim_os_sl   = sim_os_sl;                                           % exponential decay offset and slope parameter values for agent simulation (BB-H-E)
beh_inps.beta_int    = [beta_min beta_max];                                 % softmax beta interval endpoints (all agents except BF-0)     
beh_inps.beta_init   = beta_init;                                           % softmax beta initial point (all agents except BF-0)
beh_inps.lambda_int  = [lambda_min lambda_max];                             % linear convex combination lambda interval endpoints (BB-H-C)   
beh_inps.lambda_init = lambda_init;                                         % linear convex combination lambda initial point (BB-H-C)
beh_inps.os_int      = [os_min os_max];                                     % exponential decay offset interval endpoints (BB-H-E)
beh_inps.os_init     = os_init;                                             % exponential decay offset initial point (BB-H-E)  
beh_inps.sl_int      = [sl_min sl_max];                                     % exponential decay slope interval endpoints (BB-H-E)
beh_inps.sl_init     = sl_init;                                             % exponential decay slope initial point (BB-H-E)
beh_inps.plot_ll     = plot_ll;                                             % flag: plot log likelihood function 
beh_inps.stoi        = stoi;                                                % flag: evaluate agent only on task with target locations of interest
beh_inps.toi         = toi;                                                 % target locations of interest

% analyze participant and agent behavior
th_beh_analysis(beh_inps); 

% -------------------------------------------------------------------------
%                       (3) fMRI data preprocessing 
% -------------------------------------------------------------------------

% specify input structure
fmri_pre_inps           = [];                                               % initialize input structure
fmri_pre_inps.data_dir  = data_dir;                                         % group-level (=main) data directory
fmri_pre_inps.func_dir  = func_dir;                                         % functional data subdirectory: participant-level
fmri_pre_inps.pre_dir   = pre_dir;                                          % fmri preprocessing subdirectory: participant-level
fmri_pre_inps.subs      = subs;                                             % BIDS style participant IDs
fmri_pre_inps.runs      = runs;                                             % participant run numbers

% preprocess fMRI data
th_fmri_prepro(fmri_pre_inps);

% -------------------------------------------------------------------------
%                          (4) fMRI data analysis 
% -------------------------------------------------------------------------

% specify input structure
fmri_inps                 = [];                                             % initialize input structure                                    
fmri_inps.data_dir        = data_dir;                                       % group-level (=main) data directory
fmri_inps.func_dir        = func_dir;                                       % functional data subdirectory: participant-level 
fmri_inps.glm_dir_l1      = glm_dir_l1;                                     % fMRI results subdirectory: participant-level
fmri_inps.glm_dir_l2      = glm_dir_l2;                                     % fMRI results subdirectory: group-level
fmri_inps.subs            = subs;                                           % BIDS style participant IDs
fmri_inps.runs            = runs;                                           % participant run numbers
fmri_inps.t_in_run        = t_in_run;                                       % number of tasks in run
fmri_inps.max_attempt     = max_attempt;                                    % maximum number of attempts on task
fmri_inps.pomdp           = pomdp;                                          % task pomdp formulation   
fmri_inps.lumin           = [];                                             % initialize luminance input structure
fmri_inps.lumin.stim_dir  = fullfile(data_dir, stim_dir);                   % full path: visual stimuli subdirectory
fmri_inps.lumin.grid_ims  = grid_fname;                                     % visual stimuli filenames: grid cell position
fmri_inps.lumin.treas_im  = treas_fname;                                    % visual stimuli filename: treasure
fmri_inps.lumin.obs_ims   = obs_fname;                                      % visual stimuli filenames: observation bars
fmri_inps.lumin.arrow_ims = arrow_fname;                                    % visual stimuli filenames: decision prompt arrows

% analyze fMRI data
th_fmri_analysis(fmri_inps)

end

