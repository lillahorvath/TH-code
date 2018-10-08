function lleval = th_data_eval(lleval)

% This function evaluates the log likelihood for a given participant's data
% set for a number of agent models of interest. It additionally computes
% rho^2, the log likelihood based pseudo-r^2 statistics for each optimizing 
% agent (all except BF-0).
%
% Inputs
%       lleval:        structure with fields
%        .pomdp        pomdp formulation
%        .y            observable data
%        .n_run:       number of runs
%        .n_task       number of tasks per run
%        .aoi          agents of interest
%        .beta_int:    beta interval endpoints 
%        .beta_init:   beta initial point 
%        .lambda_int:  lambda interval endpoints 
%        .lambda_init: lambda initial point 
%        .os_int:      decay offset interval endpoints 
%        .os_init      decay offset initial point 
%        .sl_int:      decay slope interval endpoints 
%        .sl_init      decay slope initial point 
%        .br_sd:       subject-specific behavioral results (ll figures) directory
%        .pl:          flag:plot log likelihood function
%
% Outputs
%       lleval:        input structure with additional fields
%        .ll:          agent-specific choice data maximum log likelihoods
%        .beta_ml:     beta ML estimates
%        .lambda_ml:   lambda ML estimate
%        .os_sl_ml:    decay offset and slope ML estimates 
%        .rho_square:  log likelihood-based pseudo-r^2 statistics
%   
% Copyright (C) Dirk Ostwald, Lilla Horvath
% -------------------------------------------------------------------------

% initialization
% -------------------------------------------------------------------------

% unpack input structure
pomdp       = lleval.pomdp;
y           = lleval.y;
n_run       = lleval.n_run;
n_task      = lleval.n_task;
aoi         = lleval.aoi;
sdir        = lleval.br_sd;
beta_int    = lleval.beta_int;
beta_init   = lleval.beta_init;
lambda_int  = lleval.lambda_int;
lambda_init = lleval.lambda_init;
os_int      = lleval.os_int;
os_init     = lleval.os_init; 
sl_int      = lleval.sl_int;
sl_init     = lleval.sl_init;
plot_ll     = lleval.pl;

% initialize arrays
ll_all     = NaN(numel(aoi),1);                                             % log likelihood array
b_ml_all   = NaN(numel(aoi),1);                                             % maximum likelihood estimates - beta (softmax function, all except BF-0)
l_ml       = NaN;                                                           % maximum likelihood estimate - lambda (constant linear combination function, BB-H-C only)
os_ml      = NaN;                                                           % maximum likelihood estimate - offset (exponential decay function, BB-H-E only)
sl_ml      = NaN;                                                           % maximum likelihood estimate - slope (exponential decay function, BB-H-C only)
rho_square = NaN(numel(aoi),1);                                             % pseudo-r^2 statistics

% evaluate agents
% -------------------------------------------------------------------------

% cycle over agents
for a = 1:numel(aoi)
    
    % specify agent log likelihood maximization input structure
    a_fit             = [];                                                 % initialize structure
    a_fit.pomdp       = pomdp;                                              % pomdp formulation
    a_fit.y           = y;                                                  % observable data
    a_fit.nruns       = n_run;                                              % number of runs
    a_fit.ntasks      = n_task;                                             % number of tasks per run
    a_fit.a_class     = aoi{a};                                             % agent class
    a_fit.beta_int    = beta_int;                                           % beta optimization interval
    a_fit.beta_init   = beta_init;                                          % beta optimization initial value
    a_fit.lambda_int  = lambda_int;                                         % lambda optimization interval
    a_fit.lambda_init = lambda_init;                                        % lambda optimization initial value
    a_fit.os_int      = os_int;                                             % decay offset optimization interval
    a_fit.os_init     = os_init;                                            % decay offset optimization initial value
    a_fit.sl_int      = sl_int;                                             % decay slope optimization interval
    a_fit.sl_init     = sl_init;                                            % decay slope optimization initial value 
    a_fit.plot_ll     = plot_ll;                                            % log likelihood function plot flag
    a_fit.sdir        = sdir;                                               % participant behavioral results (ll figures) directory 
    
    % maximize agent log likelihood 
    a_fit = th_agent_fit(a_fit);

    % record ML estimates
    if isequal(a_fit.a_class, @th_agentobj_BF_0) 
        b_ml_all(a) = NaN;                                                  % no ML beta estimate (BF-0)
    else
        b_ml_all(a) = a_fit.beta_ml;                                        % Ml beta estimate (all agents except BF-0)
        if isequal(a_fit.a_class, @th_agentobj_BB_H_C)  
            l_ml = a_fit.lambda_ml;                                         % Ml lambda estimate (BB-H-C)
        elseif isequal(a_fit.a_class, @th_agentobj_BB_H_E)
            os_ml = a_fit.os_ml;                                            % ML offset estimate (BB-H-E)
            sl_ml = a_fit.sl_ml;                                            % ML slope estimate (BB-H-E)
        end
    end
    
    % record agent maximum log likelihood
    ll_all(a) = a_fit.ll;
    
    % compute rho^2 for each optimizing agent
    if a > 1
        rho_square(a) = 1-(ll_all(a)/ll_all(1));
    end
    
end

% set the maximum likelihood estimates, log likelihood array and pseudo-r^2
% statistics
lleval.beta_ml    = b_ml_all;
lleval.lambda_ml  = l_ml; 
lleval.os_sl_ml   = [os_ml, sl_ml];
lleval.ll         = ll_all;
lleval.rho_square = rho_square;

end
 


