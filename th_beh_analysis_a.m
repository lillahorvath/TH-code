function th_beh_analysis_a(beh_a_inps)

% This function evaluates the behavior of the agent models on the treasure
% hunt task. The analysis includes the following steps:
%
% (1) agent model simulations - standard task configurations
% (2) agent model recovery
% (3) agent behavior visualization
%
% Inputs
%       beh_inps_a:    structure with fields
%        .data_dir:    group-level data directory
%        .br_dir:      behavioral results subdirectory
%        .a_brn:       name of the behavioral results file - agents
%        .max_attempt: attempt limit
%        .pomdp:       task PoMDP formulation
%        .agents:      agents of interest
%        .sim_rep:     number of simulation repeats
%        .par_beta:    softmax beta parameter (all agents except BF-0)
%        .par_lambda:  linear convex combination lambda parameter (BB-H-C)
%        .par_os_sl:   exponential decay offset and slope parameter values (BB-H-E)
%        .beta_int:    beta interval endpoints - ll optimization
%        .beta_init:   beta initial point - ll optimization
%        .lambda_int:  lambda interval endpoints - ll optimization
%        .lambda_init: lambda initial point - ll optimization
%        .os_int:      decay offset interval endpoints - ll optimization
%        .os_init      decay offset initial point - ll optimization
%        .sl_int:      decay slope interval endpoints - ll optimization
%        .sl_init      decay slope initial point - ll optimization
%        .plot_ll:     flag: plot log likelihood function 
%        .toi:         target locations of interest
%        .stoi:        flag: evaluate agent only on task with target locations of interest 
%
% Outputs
%       none, writes results to disk
%
% Copyright (C) Dirk Ostwald, Lilla Horvath
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%                             initialization 
% -------------------------------------------------------------------------

% unpack input structure
data_dir    = beh_a_inps.data_dir;                                         
br_dir      = beh_a_inps.br_dir;  
a_brn       = beh_a_inps.a_brn;
natt        = beh_a_inps.max_attempt;                                       
pomdp       = beh_a_inps.pomdp;                                             
agents      = beh_a_inps.agents;                                            
nrep        = beh_a_inps.sim_rep;                                           
par_beta    = beh_a_inps.par_beta;                                          
par_lambda  = beh_a_inps.par_lambda;                                        
par_os_sl   = beh_a_inps.par_os_sl;                                         
beta_int    = beh_a_inps.beta_int;
beta_init   = beh_a_inps.beta_init;
lambda_int  = beh_a_inps.lambda_int;
lambda_init = beh_a_inps.lambda_init;
os_int      = beh_a_inps.os_int;
os_init     = beh_a_inps.os_init; 
sl_int      = beh_a_inps.sl_int;
sl_init     = beh_a_inps.sl_init;
plot_ll     = beh_a_inps.plot_ll;                                           
toi         = beh_a_inps.toi;                                               
sim_toi     = beh_a_inps.stoi;                                              
nags        = numel(agents);                                                % number of agent models
L1          = pomdp.L1;                                                     % L1 distances
nn          = pomdp.nn;                                                     % number of nodes

% create target location array
tgtloc = [];
for t1 = 2:nn
    for t2 = t1+1:nn
        tgtloc = [tgtloc; [t1 t2]];
    end
end
ntgt = size(tgtloc,1);                                                      % number of target configurations

% initialize arrays
y     = cell(ntgt,nrep,nags);                                               % observable data
x     = cell(ntgt,nrep,nags);                                               % unobservable data
B     = cell(ntgt,nrep,nags);                                               % belief state probability maps
M_sol = NaN(ntgt,nrep,nags);                                                % agent performance array (solved tasks)
M_rec = NaN(ntgt,nrep,nags,nags);                                           % agent model recovery array (data log likelihoods)
n_dp  = NaN(ntgt,nrep,nags,nags);                                           % number of valid data points array

% if evaluate agents only on task with target locations of interest 
if sim_toi                                                            
    tgtloc = [];
    tgtloc = toi;
end

% -------------------------------------------------------------------------
%                          analyze agent behavior   
% -------------------------------------------------------------------------
    
% cycle over target locations
for t = 1:size(tgtloc,1) 

    % simulation parameters
    tgt = tgtloc(t,:);                                                      % treasure locations
    ops = min(L1(1,tgt(1)),L1(1,tgt(2))) + L1(tgt(1),tgt(2));               % optimal step number
    avs = ones(1,natt)*ops;                                                 % available step number (= optimal step size (because standard task configurations))
       
    %rng('shuffle')
    
    % cycle over simulation repeats
    for r = 1:nrep    

        % inform user
        fprintf('Agent simulation and recovery - treasure configuration #%d of %d, simulation %d of %d\n', t, size(tgtloc,1), r, nrep)
      
        % cycle over agents
        for s = 1:nags

            % ******************* (1) agent simulation ********************
            % simulate agent behavior on standard task configurations
            
            % specify input structure
            a_sim            = [];                                          % initialize 
            a_sim.pomdp      = pomdp;                                       % pomdp formulation
            a_sim.tgt        = tgt;                                         % target location
            a_sim.avs        = avs;                                         % available steps
            a_sim.natt       = natt;                                        % maximal number of attempts
            a_sim.par_beta   = par_beta;                                    % beta 
            a_sim.par_lambda = par_lambda;                                  % lambda (BB-H-C)
            a_sim.par_os     = par_os_sl(1);                                % offset (BB-H-E)
            a_sim.par_sl     = par_os_sl(2);                                % slope (BB-H-E)
            a_sim.a_class    = agents{s};                                   % simulation agent class
            
            % simulate agent
            a_sim = th_agent_sim(a_sim);                          

            % record simulation
            y{t}{r}{s}   = a_sim.y;                                         % observable data
            x{t}{r}{s}   = a_sim.x;                                         % unobservable data
            B{t}{r}{s}   = a_sim.B;                                         % belief state
            M_sol(t,r,s) = a_sim.sol;                                       % agent performance indicator

            % ******************* (2) agent evaluation ********************
            % model recovery analysis: evaluate agents in light of the s-th
            % agent's simulated data
            
            % cycle over agents
            for e = 1:nags

                % specify input structure
                a_fit             = [];                                     % initialize 
                a_fit.pomdp       = pomdp;                                  % pomdp formulation
                a_fit.y{1}{1}     = y{t}{r}{s};                             % observable data
                a_fit.nruns       = 1;                                      % number of runs
                a_fit.ntasks      = 1;                                      % number of tasks per run
                a_fit.a_class     = agents{e};                              % agent class
                a_fit.beta_int    = beta_int;                               % beta optimization interval
                a_fit.beta_init   = beta_init;                              % beta optimization initial value
                a_fit.lambda_int  = lambda_int;                             % lambda optimization interval
                a_fit.lambda_init = lambda_init;                            % lambda optimization initial value
                a_fit.os_int      = os_int;                                 % decay offset optimization interval
                a_fit.os_init     = os_init;                                % decay offset optimization initial value
                a_fit.sl_int      = sl_int;                                 % decay slope optimization interval
                a_fit.sl_init     = sl_init;                                % decay slope optimization initial value 
                a_fit.sdir        = [];                                     % behavioral results directory
                a_fit.plot_ll     = plot_ll;                                % log likelihood function plot flag
                
                % maximize agent log likelihood 
                a_fit = th_agent_fit(a_fit);

                % record agent maximum log likelihood and number of data
                % points (identical across agents)
                M_rec(t,r,s,e) = a_fit.ll;
                n_dp(t,r,s,e)  = a_fit.n_dp;

            end

        end
        
        % *************** (3) agent behavior visualization ****************
        % visualize agent behavior on the task with targets of interest
        
        if r == 1 && isequal(tgt,toi)
            
            % invoke Yair Altman's export_fig toolbox (usage: save figures with minimal white space)
            addpath(genpath([pwd '\export_fig']))
            
            % cycle over agents
            for s = 1:nags
               
               % specify input structure
               taskvis       = [];                                          % initialize
               taskvis.pomdp = pomdp;                                       % pomdp formulation
               taskvis.y     = y{t}{1}{s};                                  % observable data
               taskvis.x     = x{t}{1}{s};                                  % unobservable data
               taskvis.B     = B{t}{1}{s};                                  % belief state
              
               % if a hybrid agent is evaluated
               if s == 8 || s == 9                                          
                   
                   % visualize behavior attempt-wise
                   [h1, h2, h3] = th_taskvis_BB_H(taskvis);
                   
                   % save figures 
                   if ~isempty(h1)
                        fig_name = fullfile(data_dir, br_dir, [func2str(agents{s}), '_att1.pdf']);
                        export_fig(h1, fig_name);
                    end
                    
                    if ~isempty(h2)
                        fig_name = fullfile(data_dir, br_dir, [func2str(agents{s}), '_att2.pdf']);
                        export_fig(h2, fig_name);
                    end
                    
                    if ~isempty(h3)
                        fig_name = fullfile(data_dir, br_dir, [func2str(agents{s}), '_att3.pdf']);
                        export_fig(h3, fig_name);
                    end
                    
                    close all
                    clearvars 'h1' 'h2' 'h3' 
                   
               % for all other agents
               else
                   
                    % visualize behavior task-wise
                    h = th_taskvis(taskvis); 
                    
                    % save figure
                    fig_name = fullfile(data_dir, br_dir, [func2str(agents{s}), '.pdf']);
                    export_fig(h, fig_name);
                    close
                    clearvars 'h' 
               end
            end
        end
    end

    % *************************** save results ****************************
    
    sim_pars = [par_beta par_lambda par_os_sl];
    ev_pars  = [beta_int beta_init lambda_int lambda_init os_int os_init sl_int sl_init];
    
    save(fullfile(data_dir, br_dir, a_brn), 'pomdp', 'agents', 'tgtloc', 'sim_pars', 'ev_pars', 'y', 'x', 'B', 'M_sol', 'M_rec', 'n_dp')
end

% ***************************** save results ******************************

save(fullfile(data_dir, br_dir, a_brn), 'pomdp', 'agents', 'tgtloc', 'sim_pars', 'ev_pars', 'y', 'x', 'B', 'M_sol', 'M_rec', 'n_dp')

end
 
    
 
