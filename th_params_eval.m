function params = th_params_eval(params)

% This function evaluates a set of parametric regressors for the GLM-based
% analysis of the treasure hunt fMRI dataset.
%
% Inputs
%       params:         structure with fields
%        .n_run:        number of runs
%        .n_task:       number of tasks per run
%        .pomdp:        task pomdp formulation
%        .aoi:          agent of interest (behaviorally most plausible agent BB-R-25)
%        .par_beta:     group median of the ML softmax parameter (beta) estimate
%        .y:            participant observable dataset 
%
% Outputs
%       params:         input structure with additional fields
%        .KL:           trial-by-trial KL divergence
%        .Q:            trial-by-trial chosen action valence
%        .NAC:          trial-by-trial normalized action entropy
%
% Copyright (C) Dirk Ostwald, Lilla Horvath
% -------------------------------------------------------------------------

% initialization
% -------------------------------------------------------------------------

% unpack input
nruns   = params.n_run;
ntasks  = params.n_task;
pomdp   = params.pomdp;
aoi     = params.aoi;
par     = params.par_beta; 
y       = params.y;

% initialize parametric regressor arrays
BS     = cell(size(y));
KL     = cell(size(y));
NAC    = cell(size(y));
Q      = cell(size(y));

% cycle over runs
for run = 1:nruns

    ngc  = pomdp.nn;                                                        % number of grid cells 
    rhti = [];                                                              % run-specific experimental hierarchy trial counter

    % cycle over tasks 
    for task = 1:ntasks

        % evaluate belief state, chosen action valences and normalized
        % action entropy
        % -----------------------------------------------------------------
        
        % specify input structure
        a_eva         = [];                                               
        a_eva.pomdp   = pomdp;                                              % pomdp formulation
        a_eva.y       = y{run}{task};                                       % observable data
        a_eva.a_class = aoi;                                                % agent class
        a_eva.par     = par;                                                % group median of the ML softmax parameter (beta) estimate
        
        % evaluate agent
        a_eva = th_agent_eva(a_eva);                                
        
        % assign values
        B_task   = [];
        Q_task   = [];
        NAC_task = [];
        B_task   = a_eva.B_evo;
        Q_task   = a_eva.Q_evo;
        NAC_task = a_eva.NAC_evo;
 
        % evaluate belief state KL divergences
        % -----------------------------------------------------------------
        
        % number of attempts on the current task (1 minimum to 3 maximum)
        natt = size(B_task,1);    

        % cycle over attempts
        for att = 1:natt

            % extract belief states of the current attempt
            B_task_att = B_task{att};
            
            % determine number of trials, including last position presentation - the first belief state preceeds the first trial
            ntrials = numel(B_task_att)-1;
            
            % BS evolution initialization
            KL_trials = NaN(ntrials,1);

            % cycle over trials and evaluate belief state KL(posterior||prior)
            % [sequential zero allocations require this order, KL(P||Z) is defined
            % only if Z(i) implies P(i) for all i ("absolute continuity")]
            for trial = 1:ntrials

                % concatentate the experimental hierarchy trial index
                rhti = [rhti [task att trial]'];

                % initialize KL divergence
                KL_trials(trial) = 0;

                % cycle over grid cells
                for gc = 1:ngc

                    % exclude zero probability grid cells
                    if B_task_att{trial}(gc) ~= 0 && B_task_att{trial+1}(gc) ~= 0
                        KL_trials(trial) = KL_trials(trial) + B_task_att{trial+1}(gc)*log(B_task_att{trial+1}(gc)/B_task_att{trial}(gc));
                    end
                end

            end
            
            % save value evolution
            KL{run}{task}{att}  = KL_trials;
            BS{run}{task}{att}  = B_task_att;
            Q{run}{task}{att}   = Q_task{att};
            NAC{run}{task}{att} = NAC_task{att};
            
        end
    end
end

% assign to output structure
params.KL   = KL;
params.Q    = Q;
params.NAC  = NAC;

end



