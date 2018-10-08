function th_beh_analysis_p(beh_p_inps)

% This function evaluates the participant behavioral data acquired on the 
% treasure hunt task. The analysis includes the following steps:
%
% (1) data conversion from BIDS events long format to wide format
% (2) descriptive statistical analysis
% (3) agent model simulations - participant task configurations
% (4) agent model evaluation 
% (5) Bayesian model selection
%
% Inputs
%       beh_inps_p:    structure with fields
%        .data_dir:    group-level data directory
%        .func_dir:    functional data subdirectory: participant-level
%        .br_dir:      behavioral results subdirectory: group-level
%        .p_brn:       name of the behavioral results file - participants
%        .subs:        BIDS style subject IDs
%        .runs:        subject run numbers
%        .max_run:     maximum number of runs
%        .t_in_run:    number of tasks in run
%        .max_attempt: maximum number of attempts on task
%        .max_opt_dis: maximum length of shortest path (=optimal step size)
%        .sln_set:     set size of the step limit noise
%        .dim:         grid-world dimensionality
%        .pomdp:       task pomdp formulation
%        .agents:      agents of interest
%        .nrep:        number of simulation repeats
%        .par_beta:    beta parameter value for agent simulations
%        .par_lambda:  lambda parameter value for agent simulation
%        .par_os_sl:   offset and slope parameter values for agent simulation
%        .beta_int:    beta interval endpoints - ll optimization
%        .beta_init:   beta initial point - ll optimization
%        .lambda_int:  lambda interval endpoints - ll optimization
%        .lambda_init: lambda initial point - ll optimization
%        .os_int:      decay offset interval endpoints - ll optimization
%        .os_init      decay offset initial point - ll optimization
%        .sl_int:      decay slope interval endpoints - ll optimization
%        .sl_init      decay slope initial point - ll optimization
%        .plot_ll:     flag: plot log likelihood function
%
% Outputs
%       none, saves data directly to disk
%
% Copyright (C) Dirk Ostwald, Lilla Horvath
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%                             initialization 
% -------------------------------------------------------------------------

% unpack input structure
data_dir    = beh_p_inps.data_dir;
func_dir    = beh_p_inps.func_dir; 
br_dir      = beh_p_inps.br_dir; 
p_brn       = beh_p_inps.p_brn;
subs        = beh_p_inps.subs;        
runs        = beh_p_inps.runs;        
max_run     = beh_p_inps.max_run;     
t_in_run    = beh_p_inps.t_in_run;    
max_attempt = beh_p_inps.max_attempt; 
max_opt_dis = beh_p_inps.max_opt_dis;
sln_set     = beh_p_inps.sln_set;     
dim         = beh_p_inps.dim;         
pomdp       = beh_p_inps.pomdp;  
aoi         = beh_p_inps.agents;                                                            
nrep        = beh_p_inps.sim_rep;                                        
par_beta    = beh_p_inps.par_beta;
par_lambda  = beh_p_inps.par_lambda;                          
par_os_sl   = beh_p_inps.par_os_sl;
beta_int    = beh_p_inps.beta_int;
beta_init   = beh_p_inps.beta_init;
lambda_int  = beh_p_inps.lambda_int;
lambda_init = beh_p_inps.lambda_init;
os_int      = beh_p_inps.os_int;
os_init     = beh_p_inps.os_init; 
sl_int      = beh_p_inps.sl_int;
sl_init     = beh_p_inps.sl_init;
plot_ll     = beh_p_inps.plot_ll;

% group data arrays
agents_all     = NaN(numel(subs), numel(aoi));                              % agent performance
models_all     = NaN(numel(subs), numel(aoi));                              % agent model log likelihoods
beta_ml_all    = NaN(numel(subs), numel(aoi));                              % beta ml estimates (softmax)
lambda_ml_all  = NaN(numel(subs), 1);                                       % ml estimate of the @th_agentobj_BB_H_C linear combination parameter
os_sl_ml_all   = NaN(numel(subs), 2);                                       % ml estimates of the @th_agentobj_BB_H_E offset and   parameter
rho_square_all = NaN(numel(subs), numel(aoi));                              % pseudo-r^2 statistics
qa_1           = NaN(numel(subs), dim*dim, dim*dim);                        % number of tasks per treasure location configurations 
qa_2           = NaN(numel(subs), max_opt_dis+1);                           % number of tasks per optimal step size
qa_3           = NaN(numel(subs), sln_set);                                 % number of attempts per step limit noise
bp_1           = NaN(numel(subs), 3);                                       % experiment task statistics
bp_2           = NaN(numel(subs), max_run);                                 % run performance (percentage, solved out of solvable)
bp_3a          = NaN(numel(subs), max_opt_dis+1);                           % performance per opt step size (percentage, solved out of solvable)
bp_3b          = NaN(numel(subs), dim*dim, dim*dim);                        % performance per treasure loc configs (percentage, solved out of solvable)
brt_1          = NaN(numel(subs), 1);                                       % median RT
brt_2          = NaN(numel(subs), max_run);                                 % median RT per run
brt_3          = NaN(numel(subs), max_opt_dis+1);                           % median RT per opt step size
ap_1           = NaN(numel(subs), max_attempt);                             % performance per attempt (percentage, solved in attempt out of all solvable)
ap_2           = NaN(numel(subs), max_attempt);                             % performance per attempt (percentage, solved in attempt out of all solved)
art            = NaN(numel(subs), max_attempt);                             % median RT per attempt 
dec_nums       = NaN(numel(subs), 6);                                       % number of decisions: total num, missing or not valid, towards a dark bar, towards a light bar, towards a db in the presence of a lb, towards a lb in the presence of a db 
bic_all        = NaN(numel(subs), numel(aoi));                              % agent model BIC scores
av_pp          = NaN(numel(aoi), 1);                                        % hypothesis-averaged posterior model probabilities

% -------------------------------------------------------------------------
%                    analyze participant behavioral data   
% -------------------------------------------------------------------------

% cycle over participants
for sj = 1:numel(subs)
    
    % **************************** initialize *****************************
    
    % inform user
    fprintf(['Analyzing ' subs{sj}, '\n'])
    
    % ************************* (1) convert data **************************
    % convert data from BIDS long format to wide format: create cell
    % arrays y and x containing the run-task-attempt observable and not 
    % observable data 
    
    % specify input structure 
    data_yx         = [];                                                   % initialize structure
    data_yx.sub_dir = fullfile(data_dir, subs{sj}, func_dir);               % full path: participant-level functional data subdirectory
    data_yx.sub     = subs{sj};                                             % BIDS style subject ID
    data_yx.sub_run = runs{sj};                                             % subject run 
    data_yx.n_task  = t_in_run;                                             % number of tasks in run
    data_yx.a_max   = max_attempt;                                          % maximum number of attempts
    
    % convert data 
    data_yx = th_beh_data_yx(data_yx);               
 
    % ******************** (2) descriptive statistics *********************
    % compute summary statistics of the participant's data 
    
    % specify input structure
    stats         = [];                                                     % initialize structure 
    stats.sub_run = runs{sj};                                               % subject run
    stats.pomdp   = pomdp;                                                  % pomdp formulation  
    stats.n_task  = t_in_run;                                               % number of tasks in run
    stats.a_max   = max_attempt;                                            % max num of attempts
    stats.sp_max  = max_opt_dis;                                            % maximum length of shortest path (=optimal step size)
    stats.y       = data_yx.y;                                              % observable data
    stats.x       = data_yx.x;                                              % unobservable data

    % obtain task statistics
    stats = th_data_stats(stats);
    
    % *********************** (3) agent simulation ************************
    % simulate agent behavior on the participant's task configurations
    
    % specify input struture
    tasks           = [];                                                   % initialize structure
    tasks.agents    = aoi;                                                  % agents of interest 
    tasks.taskall   = stats.tasks;                                          % task configurations
    tasks.pomdp     = pomdp;                                                % pomdp formulation
    tasks.nrep      = nrep;                                                 % number of simulation repeats
    tasks.par_be    = par_beta;                                             % beta parameter value (all agents except BF-0)
    tasks.par_la    = par_lambda;                                           % lambda parameter value (BB-H-C)
    tasks.par_os_sl = par_os_sl;                                            % offset and slope parameter values (BB-H-E)
    
    % simulate agents on participant task configuration
    tasks = th_data_sim(tasks);
     
    % *********************** (4) agent evaluation ************************
    % evaluate agents in light of the participant's data 
    
    % specify input structure
    lleval             = [];                                                % initialize structure 
    lleval.n_run       = numel(runs{sj});                                   % number of runs
    lleval.n_task      = t_in_run;                                          % number of tasks per run
    lleval.pomdp       = pomdp;                                             % pomdp formulation
    lleval.aoi         = aoi;                                               % agents of interest
    lleval.y           = data_yx.y;                                         % observable data
    lleval.beta_int    = beta_int;                                          % beta interval for ll optimization
    lleval.beta_init   = beta_init;                                         % beta initial value for ll optimization
    lleval.lambda_int  = lambda_int;                                        % lambda interval for ll optimization
    lleval.lambda_init = lambda_init;                                       % lambda initial value for ll optimization
    lleval.os_int      = os_int;                                            % decay offset interval for ll optimization
    lleval.os_init     = os_init;                                           % decay offset initial value for ll optimization
    lleval.sl_int      = sl_int;                                            % decay slope interval for ll optimization
    lleval.sl_init     = sl_init;                                           % decay slope initial value for ll optimization
    lleval.pl          = plot_ll;                                           % flag: plot log likelihood function
    lleval.br_sd       = fullfile(data_dir, subs{sj}, [br_dir '_' subs{sj}]); % full path: participant-level behavioral results (ll figures) subdirectory

    % evaluate agent-specific data log likelihood
    lleval = th_data_eval(lleval);
    
    % *************************** allocate data ***************************
    
    qa_1(sj,:,:)                = stats.tl(:,:,1)';
    qa_2(sj,:)                  = stats.dist(:,2);
    qa_3(sj,:)                  = stats.step_lim; 
    bp_1(sj,:)                  = [stats.ntask stats.solvable stats.solved];
    bp_2(sj,1:numel(runs{sj}))  = stats.run(:,4);
    bp_3a(sj,:)                 = stats.dist(:,5);
    bp_3b(sj,:,:)               = stats.tl(:,:,4)';
    brt_1(sj)                   = stats.RT.avg;
    brt_2(sj,1:numel(runs{sj})) = stats.run_rt;
    brt_3(sj,:)                 = stats.dist_rt;
    ap_1(sj,:)                  = stats.attempt(:,1);
    ap_2(sj,:)                  = stats.attempt(:,2);
    art(sj,:)                   = stats.attempt_rt;
    dec_nums(sj,:)              = stats.dn;   
    agents_all(sj,:)            = squeeze(mean(sum(tasks.M_sol,1),2))';
    models_all(sj,:)            = lleval.ll;
    beta_ml_all(sj,:)           = lleval.beta_ml;
    lambda_ml_all(sj)           = lleval.lambda_ml;
    os_sl_ml_all(sj,:)          = lleval.os_sl_ml; 
    rho_square_all(sj,:)        = lleval.rho_square;
    
    % *************************** save results ****************************
    
    sim_pars = [par_beta par_lambda par_os_sl];
    ev_pars  = [beta_int beta_init lambda_int lambda_init os_int os_init sl_int sl_init];
    
    save(fullfile(data_dir, br_dir, p_brn), 'subs', 'pomdp', 'aoi', 'max_run', 'max_attempt', 'max_opt_dis', 'sln_set', 'sim_pars', 'ev_pars', 'dec_nums', 'bp_1', 'bp_2', 'bp_3a', 'bp_3b', 'brt_1', 'brt_2', 'brt_3', 'art', 'ap_1', 'qa_1', 'qa_2', 'qa_3', 'agents_all', 'beta_ml_all', 'lambda_ml_all', 'os_sl_ml_all', 'rho_square_all', 'bic_all', 'av_pp');
    
end

% ******************* (5) Bayesian model selection ********************
% evaluate the agent model exceedance probabilities based on the BIC scores

% specify input structure
bmselect    = [];
bmselect.ll = models_all;                                                   % agent models maximum log likelihood - all participants
bmselect.dn = dec_nums(:,1);                                                % number of valid decisions - all participants

% perform Bayesian model selection
bmselect = th_agent_select(bmselect);

% ***************************** allocate data *****************************

bic_all = bmselect.bic;
av_pp   = bmselect.pxp;

% ***************************** save results ******************************

save(fullfile(data_dir, br_dir, p_brn), 'subs', 'pomdp', 'aoi', 'max_run', 'max_attempt', 'max_opt_dis', 'sln_set', 'sim_pars', 'ev_pars', 'dec_nums', 'bp_1', 'bp_2', 'bp_3a', 'bp_3b', 'brt_1', 'brt_2', 'brt_3', 'art', 'ap_1', 'qa_1', 'qa_2', 'qa_3', 'agents_all', 'beta_ml_all', 'lambda_ml_all', 'os_sl_ml_all', 'rho_square_all', 'bic_all', 'av_pp');

end 




