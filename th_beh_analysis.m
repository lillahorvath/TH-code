function th_beh_analysis(beh_inps)

% This function evaluates the behavior of the participants and agent models 
% on the treasure hunt task. It performs the following steps:
% (1) behavioral data analysis - participant data 
% (2) behavioral data analysis - agents
% (3) results visualization
% 
% Inputs
%       beh_inps:      structure with fields
%        .data_dir:    group-level data directory
%        .func_dir:    functional data subdirectory: participant-level
%        .br_dir:      behavioral results subdirectory: group-level 
%        .p_brn:       name of the behavioral results file - participants
%        .a_brn:       name of the behavioral results file - agents
%        .subs:        BIDS style participant IDs
%        .runs:        participant run numbers
%        .max_run      maximum number of experimental runs 
%        .t_in_run:    number of tasks in run 
%        .max_attempt: maximum number of attempts on task
%        .max_opt_dis: maximum length of shortest path (= optimal step size) 
%        .sln_set:     set size of the step limit noise 
%        .dim:         grid-world dimensionality         
%        .pomdp:       task pomdp formulation
%        .agents:      agents of interest                            
%        .nrep:        number of simulation repeats
%        .sim_beta:    softmax beta parameter for agent simulations (all agents except BF-0)
%        .sim_lambda:  linear convex combination lambda parameter for agent simulation (BB-H-C)
%        .sim_os_sl:   exponential decay offset and slope parameter values for agent simulation (BB-H-E)
%        .beta_int:    softmax beta interval endpoints (all agents except BF-0)     
%        .beta_init:   softmax beta initial point (all agents except BF-0)
%        .lambda_int:  linear convex combination lambda interval endpoints (BB-H-C)   
%        .lambda_init: linear convex combination lambda initial point (BB-H-C)
%        .os_int:      exponential decay offset interval endpoints (BB-H-E)
%        .os_init      exponential decay offset initial point (BB-H-E)  
%        .sl_int:      exponential decay slope interval endpoints (BB-H-E)
%        .sl_init      exponential decay slope initial point (BB-H-E)
%        .plot_ll:     flag: plot log likelihood function 
%        .stoi:        flag: evaluate agent only on task with target locations of interest
%        .toi:         target locations of interest
%
% Outputs
%       none
%
% Copyright (C) Lilla Horvath, Dirk Ostwald
% -------------------------------------------------------------------------

% initialize
% -------------------------------------------------------------------------

% unpack input structure
data_dir    = beh_inps.data_dir;
func_dir    = beh_inps.func_dir; 
br_dir      = beh_inps.br_dir;
p_brn       = beh_inps.p_brn;
a_brn       = beh_inps.a_brn;
subs        = beh_inps.subs;
runs        = beh_inps.runs;
max_run     = beh_inps.max_run;
t_in_run    = beh_inps.t_in_run;
max_attempt = beh_inps.max_attempt;
max_opt_dis = beh_inps.max_opt_dis;
sln_set     = beh_inps.sln_set;
dim         = beh_inps.dim;
pomdp       = beh_inps.pomdp;
agents      = beh_inps.agents;                       
nrep        = beh_inps.nrep;
sim_beta    = beh_inps.sim_beta;
sim_lambda  = beh_inps.sim_lambda;
sim_os_sl   = beh_inps.sim_os_sl;
beta_int    = beh_inps.beta_int;
beta_init   = beh_inps.beta_init;
lambda_int  = beh_inps.lambda_int;
lambda_init = beh_inps.lambda_init;
os_int      = beh_inps.os_int;
os_init     = beh_inps.os_init; 
sl_int      = beh_inps.sl_int;
sl_init     = beh_inps.sl_init;
plot_ll     = beh_inps.plot_ll;
stoi        = beh_inps.stoi;
toi         = beh_inps.toi;

% create behavioral results subdirectory
if exist(fullfile(data_dir, br_dir), 'dir')                                 % delete if pre-existent and create a new one
    rmdir(fullfile(data_dir, br_dir),'s')
    mkdir(fullfile(data_dir, br_dir))
else
    mkdir(fullfile(data_dir, br_dir))                                       % create if non-existent
end

% analyze behavior 
% -------------------------------------------------------------------------

% --------------------------- participant data ----------------------------

% specify input structure
beh_p_inps             = [];                                                % initialize input structure
beh_p_inps.data_dir    = data_dir;                                          % group-level data directory
beh_p_inps.func_dir    = func_dir;                                          % functional data subdirectory: participant-level
beh_p_inps.br_dir      = br_dir;                                            % behavioral results subdirectory: group-level
beh_p_inps.p_brn       = p_brn;                                             % name of the behavioral results file - participants
beh_p_inps.subs        = subs;                                              % BIDS style participant IDs
beh_p_inps.runs        = runs;                                              % participant run numbers
beh_p_inps.max_run     = max_run;                                           % maximum number of experimental runs 
beh_p_inps.t_in_run    = t_in_run;                                          % number of tasks in run 
beh_p_inps.max_attempt = max_attempt;                                       % maximum number of attempts on task
beh_p_inps.max_opt_dis = max_opt_dis;                                       % maximum length of shortest path (= optimal step size) 
beh_p_inps.sln_set     = sln_set;                                           % set size of the step limit noise 
beh_p_inps.dim         = dim;                                               % grid-world dimensionality         
beh_p_inps.pomdp       = pomdp;                                             % task pomdp formulation
beh_p_inps.agents      = agents;                                            % agents of interest                            
beh_p_inps.sim_rep     = nrep;                                              % number of simulation repeats
beh_p_inps.par_beta    = sim_beta;                                          % softmax beta parameter for agent simulations (all agents except BF-0)
beh_p_inps.par_lambda  = sim_lambda;                                        % linear convex combination lambda parameter for agent simulation (BB-H-C)
beh_p_inps.par_os_sl   = sim_os_sl;                                         % exponential decay offset and slope parameter values for agent simulation (BB-H-E)
beh_p_inps.beta_int    = beta_int;                                          % softmax beta interval endpoints (all agents except BF-0)     
beh_p_inps.beta_init   = beta_init;                                         % softmax beta initial point (all agents except BF-0)
beh_p_inps.lambda_int  = lambda_int;                                        % linear convex combination lambda interval endpoints (BB-H-C)   
beh_p_inps.lambda_init = lambda_init;                                       % linear convex combination lambda initial point (BB-H-C)
beh_p_inps.os_int      = os_int;                                            % exponential decay offset interval endpoints (BB-H-E)
beh_p_inps.os_init     = os_init;                                           % exponential decay offset initial point (BB-H-E)  
beh_p_inps.sl_int      = sl_int;                                            % exponential decay slope interval endpoints (BB-H-E)
beh_p_inps.sl_init     = sl_init;                                           % exponential decay slope initial point (BB-H-E)
beh_p_inps.plot_ll     = plot_ll;                                           % flag: plot log likelihood function 

% analyze participant behavioral data
th_beh_analysis_p(beh_p_inps); 

% ----------------------------- agent models ------------------------------

% specify input structure
beh_a_inps             = [];                                                % initialize input structure
beh_a_inps.data_dir    = data_dir;                                          % group-level data directory
beh_a_inps.br_dir      = br_dir;                                            % behavioral results subdirectory: group-level 
beh_a_inps.a_brn       = a_brn;                                             % name of the behavioral results file - agents
beh_a_inps.max_attempt = max_attempt;                                       % maximum number of attempts on task        
beh_a_inps.pomdp       = pomdp;                                             % task pomdp formulation   
beh_a_inps.agents      = agents;                                            % agents of interest
beh_a_inps.sim_rep     = nrep;                                              % number of simulation repeats
beh_a_inps.par_beta    = sim_beta;                                          % softmax beta parameter for agent simulations (all agents except BF-0)
beh_a_inps.par_lambda  = sim_lambda;                                        % linear convex combination lambda parameter for agent simulation (BB-H-C)
beh_a_inps.par_os_sl   = sim_os_sl;                                         % exponential decay offset and slope parameter values for agent simulation (BB-H-E)
beh_a_inps.beta_int    = beta_int;                                          % softmax beta interval endpoints (all agents except BF-0)     
beh_a_inps.beta_init   = beta_init;                                         % softmax beta initial point (all agents except BF-0)
beh_a_inps.lambda_int  = lambda_int;                                        % linear convex combination lambda interval endpoints (BB-H-C)   
beh_a_inps.lambda_init = lambda_init;                                       % linear convex combination lambda initial point (BB-H-C)
beh_a_inps.os_int      = os_int;                                            % exponential decay offset interval endpoints (BB-H-E)
beh_a_inps.os_init     = os_init;                                           % exponential decay offset initial point (BB-H-E)  
beh_a_inps.sl_int      = sl_int;                                            % exponential decay slope interval endpoints (BB-H-E)
beh_a_inps.sl_init     = sl_init;                                           % exponential decay slope initial point (BB-H-E)
beh_a_inps.plot_ll     = plot_ll;                                           % flag: plot log likelihood function 
beh_a_inps.toi         = toi;                                               % target locations of interest
beh_a_inps.stoi        = stoi;                                              % flag: evaluate agent only on task with target locations of interest 

% analyze agent behavior
th_beh_analysis_a(beh_a_inps);

% --------------------------- visualize results ---------------------------

% invoke Yair Altman's export_fig toolbox (usage: save figures with minimal white space)
addpath(genpath([pwd '\export_fig']))

% specify figures input structure
fig_inps          = [];                                                     % initialize
fig_inps.data_dir = data_dir;                                               % group-level data directory
fig_inps.br_dir   = br_dir;                                                 % behavioral results subdirectory: group-level
fig_inps.p_brn    = p_brn;                                                  % name of the behavioral results file - participants
fig_inps.a_brn    = a_brn;                                                  % name of the behavioral results file - agents

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% TO DO: THESE SCRIPTS ALL STORE THE FIGS WITH THE EXPORT FIG TOOLBOX AND
% WITH SIZES ADJUSTED FOR THE PAPER. THIS MIGHT NEED FURTHER ADJUSTMENTS
% WHEN FINALIZING THE PAPER
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% generate Figure 3 showing the main behavioral results
th_fig_3(fig_inps)

% generate Supplementary Figure 10 showing the basic properties of the tasks
% presented to the participants
th_sfig_10(fig_inps)

% generate Supplementary Figure 11 showing supplementary participant
% behavioral results
th_sfig_11(fig_inps)

% generate Supplementary Figure 12 showing the RT data
th_sfig_12(fig_inps)

% generate Supplementary Figure 13 showing the performance of the agents
% for each treasure location combinations
th_sfig_13(fig_inps)

% generate Supplementary Figure 14 showing the model recovery results
th_sfig_14(fig_inps)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% END OF THE TO DO
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end



