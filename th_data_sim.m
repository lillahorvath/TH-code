function tasks = th_data_sim(tasks)

% This function simulates the behavior of the agent models on the task 
% configurations presented to the participants.
%
% Inputs
%       tasks:        input structure with fields
%        .tasksall:   task configurations of a single participant
%        .pomdp:      pomdp formulation
%        .agents:     cell array of function handles with agents of interest
%        .nrep:       number of simulation repeats
%        .beta_sim:   value of beta parameter for simulation (all agents except BF-0)
%        .lambda_sim: value of lambda parameter for simulation (BB-H-C)
%        .os_sl_sim:  value of offset and slope parameters for simulation (BB-H-E)
%
% Outputs
%       tasks:        input structure with additional fields
%        .M_sol:      total attempts x repeats x agents success array
%
% Copyright (C) Dirk Ostwald, Lilla Horvath
% -------------------------------------------------------------------------

% initialize
% -------------------------------------------------------------------------

% unpack input structure
taskall    = tasks.taskall;                                            
pomdp      = tasks.pomdp;                                              
agents     = tasks.agents;                                             
nrep       = tasks.nrep;                                               
beta_sim   = tasks.par_be;                                             
lambda_sim = tasks.par_la;                                             
os_sl_sim  = tasks.par_os_sl;                                          
    
% simulation parameters
ntsk = numel(taskall);                                                      % number of task configurations
nags = numel(agents);                                                       % number of agent models

% simulation initialization
y     = cell(ntsk,nrep,nags);                                               % observable data
x     = cell(ntsk,nrep,nags);                                               % unobservable data
B     = cell(ntsk,nrep,nags);                                               % belief state probability maps
M_sol = NaN(ntsk,nrep,nags);                                                % agent model performance array (solved tasks)

% simulate agent behavior
% -------------------------------------------------------------------------

% cycle over task configurations
for t = 1:ntsk

    % cycle over agents
    for s = 1:nags
        
        % specify input structure
        a_sim            = [];                                              % initialize simulation structure
        a_sim.natt       = size(taskall{t},2);                              % number of participant's attempts
        a_sim.pomdp      = pomdp;                                           % pomdp formulation
        a_sim.tgt        = taskall{t}(1:2,1)';                              % treasure locations, constant over the attempts of a task
        a_sim.avs        = taskall{t}(4,:);                                 % available steps per attempt
        a_sim.a_class    = agents{s};                                       % agent specific decision rule
        a_sim.par_beta   = beta_sim;                                        % softmax beta parameter
        a_sim.par_lambda = lambda_sim;                                      % lambda parameter (BB-H-C)
        a_sim.par_os     = os_sl_sim(1);                                    % offset parameter (BB-H-E)
        a_sim.par_sl     = os_sl_sim(2);                                    % slope parameter (BB-H-E)
      
        % cycle over simulation repeats
        for r = 1:nrep            
   
            % simulate agent
            a_sim = th_agent_sim(a_sim);
 
            % record simulation
            y{t}{r}{s}   = a_sim.y;                                         % observable data
            x{t}{r}{s}   = a_sim.x;                                         % unobservable data
            B{t}{r}{s}   = a_sim.B;                                         % belief state
            M_sol(t,r,s) = a_sim.sol;                                       % agent performance indicator

        end 
    end
end
    
% set output structure
% -------------------------------------------------------------------------
tasks.M_sol = M_sol;                                                        % overall agent performance indicator

end