function [a_sim] = th_agent_sim(a_sim)

% This function evaluates the performance of an agent model implemented 
% using object oriented programming on the treasure hunt task. Different 
% agents belong to different agent classes. 
%
% Inputs
%       a_sim:        structure with fields
%        .pomdp:      pomdp task formulation 
%        .tgt:        1 x 2 array of treasure nodes
%        .natt:       scalar, maximum number of attempts
%        .avs:        1 x natt, attempt step limit
%        .a_class:    agent class
%        .par_beta:   softmax beta 
%        .par_lambda: BB-H-C lambda
%        .par_os:     BB-H-E offset
%        .par_sl:     BB-H-E slope
%           
% Outputs
%       a_sim:        input structure with additional fields
%        .y:          observable data
%        .x:          unobservable data
%        .B:          belief state
%        .sol:        task solved flag
%               
% Copyright (C) Dirk Ostwald, Lilla Horvath  
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%                              initialization 
% -------------------------------------------------------------------------

% unpack input structure
nn         = a_sim.pomdp.nn;                                                % number of grid-world nodes
L1         = a_sim.pomdp.L1;                                                % L1 distance map
A          = a_sim.pomdp.A;                                                 % set of actions
O          = a_sim.pomdp.O;                                                 % observation likelihood function
natt       = a_sim.natt; 
tgt        = a_sim.tgt; 
avs        = a_sim.avs; 
par_beta   = a_sim.par_beta; 
par_lambda = a_sim.par_lambda;
par_os     = a_sim.par_os;
par_sl     = a_sim.par_sl;
a_class    = a_sim.a_class; 
                                                
% initialize data array
y = cell(natt,1);                                                           % observation & action data array
x = cell(natt,1);                                                           % state and action valence array
B = cell(natt,1);                                                           % belief state array

% ************************ initialize task object *************************

% specify task initialization structure
t_init.nn  = nn;                                                            % number of gridworld nodes
t_init.L1  = L1;                                                            % L1 distances
t_init.A   = A;                                                             % set of actions
t_init.O   = O;                                                             % observation bar likelihood function
t_init.sol = 0;                                                             % initialize task solved flag
t_init.tkn = NaN;                                                           % initialize known target location
t_init.att = 0;                                                             % initialize attempt counter
t_init.t   = NaN;                                                           % initialize trial counter

% initialize task
task = th_taskobj(t_init);                                                  % task object creation

% ************************ initialize agent object ************************

% specify agent initialization structure
a_init.nn  = nn;                                                            % number of gridworld nodes
a_init.L1  = L1;                                                            % L1 distances
a_init.A   = A;                                                             % set of actions
a_init.O   = O;                                                             % observation likelihood function
a_init.be  = par_beta;                                                      % agent softmax parameter
a_init.att = 0;                                                             % initialize attempt counter
a_init.B   = [];                                                            % initialize belief state
a_init.tkn = NaN;                                                           % initialize known target location
a_init.t   = NaN;                                                           % initialize trial counter

% if agent BB-H-C is evaluated
if isequal(a_class,  @th_agentobj_BB_H_C)  
    a_init.la = par_lambda;                                                 % lambda 
% if agent BB-H-E is evaluated
elseif isequal(a_class,  @th_agentobj_BB_H_E)  
    a_init.os = par_os;                                                     % offset 
    a_init.sl = par_sl;                                                     % slope 
end

% initialize agent
agent = feval(a_class, a_init);                                             % agent object creation

% initialize task-attempt arrays
td  = cell(natt,1);                                                         % cell array with target location visits per attempts
T_a = [];                                                                   % array of known target locations from previous attempt(s)
att = 0;                                                                    % initialize attempt counter

% -------------------------------------------------------------------------
%                          simulate agent behavior  
% -------------------------------------------------------------------------

% loop over attempts
while att < natt 
    
    % ******************* initialize attempt variables ********************
    
    % update and set attempt-dependent variables
    % ---------------------------------------------------------------------
    att = att + 1;                                                          % update attempt counter
    T_a = [T_a, task.tkn];                                                  % set known target locations from previous attempt(s)
    T   = avs(att);                                                         % determine step size limit
    
    % initialize observed and unobserved random variable arrays
    % ---------------------------------------------------------------------
    s_rec = NaN(T+1,4);                                                     % state array (agent position, target found, target location)
    o_rec = NaN(T+1,4);                                                     % observation bars array
    Q_rec = NaN(T+1,size(A,2));                                             % action valence recorder
    p_rec = NaN(T+1,size(A,2));                                             % action probability recorder 
    a_rec = NaN(T+1,1);                                                     % action recorder
    
    % if agent BB-H-C is evaluated
    if isequal(a_class, @th_agentobj_BB_H_C)
        Q_rec_BB_R_1 = NaN(T+1,size(A,2));                                  % action valence recorder: BB-R-1
        Q_rec_BB_I   = NaN(T+1,size(A,2));                                  % action valence recorder: BB-I
        p_rec_BB_R_1 = NaN(T+1,size(A,2));                                  % action probability recorder: BB-R-1
        p_rec_BB_I   = NaN(T+1,size(A,2));                                  % action probability recorder: BB-I
    % if agent BB-H-E is evaluated
    elseif isequal(a_class, @th_agentobj_BB_H_E)
        Q_rec_BB_R_2 = NaN(T+1,size(A,2));                                  % action valence recorder: BB-R-2
        Q_rec_BB_I   = NaN(T+1,size(A,2));                                  % action valence recorder: BB-I
        p_rec_BB_R_2 = NaN(T+1,size(A,2));                                  % action probability recorder: BB-R-2
        p_rec_BB_I   = NaN(T+1,size(A,2));                                  % action probability recorder: BB-I
    end

    % initialize agent on attempt
    % ---------------------------------------------------------------------
    agent.t   = 0;                                                          % trial counter
    agent.T_a = T_a;                                                        % assign target found history to agent object
    agent.att = att;                                                        % attempt number
    agent.ap  = 0;                                                          % initial attempt phase
    
    % evaluate agent initial belief state
    agent = belief_update(agent);
    
    % record belief state evolution, if existent
    B_rec = [];
    
    if ~isempty(agent.B)
        B_rec      = NaN(T+2,nn);                                           % initialize attempt belief state recorder
        B_rec(1,:) = agent.B(:,end)';                                       % agent prior belief state 
    end
   
    % initialize task state on attempt
    % ---------------------------------------------------------------------
    task.t   = 1;                                                           % task trial counter
    task.att = att;                                                         % task attempt
    task.fnd = zeros(1,2);                                                  % initialize target found flag for the current attempt
    task.s_t = [1 0 tgt];                                                   % state: agent start position, attempt target found flag, target node 1, target node 2 
    
    % cycle over trials
    t = 1;                                                                  % trial index                                       
    while t <= T 
        
        % task and agent trial counter
        agent.t = t;
        task.t  = t;
        
        % ******************** simulate agent on trial ********************
        
        % task state evaluation
        % -----------------------------------------------------------------
        ft = 0;                                                             % initialize trial target found index
        if task.s_t(1) == task.s_t(3)
             task.fnd(1) = 1;                                               % first target (target at node with smaller lin. index) found
             task.tkn    = task.s_t(3);                                     % set known target location
             task.s_t(2) = 1;                                               % set attempt target found flag
             ft          = 1;                                               % turn on trial target found index
             td{att}     = [td{att}, task.s_t(3)];                          % set visited target location 
        elseif task.s_t(1) == task.s_t(4)
             task.fnd(2) = 1;                                               % second target (target at node with larger lin. index) found
             task.tkn    = task.s_t(4);                                     % set known target location
             task.s_t(2) = 2;                                               % set attempt target found flag
             ft          = 1;                                               % turn on trial target found index
             td{att}     = [td{att}, task.s_t(4)];                          % set visited target location 
        end
            
        % record state
        s_rec(t,:) = task.s_t;
        
        % if both treasures found - end the block
        if sum(task.fnd) == 2 
            task.sol = 1;
            break
        end
              
        % task observation bar sample 
        % -----------------------------------------------------------------
        task = sample_observation(task);
        
        % record the observation bars
        o_rec(t,:) = task.o_t; 
        
        % agent belief state update
        % -----------------------------------------------------------------
        agent.t   = t;                                                      % trial counter
        agent.s_t = task.s_t(1:2);                                          % observable part of the state: position and attempt target found flag
        agent.ft  = ft;                                                     % trial target found index
        agent.td  = td;                                                     % visited target locations 
        agent.o_t = task.o_t;                                               % observation bars
        agent.att = task.att;                                               % attempt on task
        agent.ap  = 1;                                                      % trial attempt phase 
        agent.tkn = task.tkn;                                               % known target location
        
        % update belief
        agent = belief_update(agent);                                
 
        % record internal agent belief state
        if ~isempty(agent.B)
            B_rec(t+1,:) = agent.B(:,end)'; 
        end
        
        % agent decision
        % -----------------------------------------------------------------
        agent = decide(agent);
       
        % record agent internal decision variables: agent action valences
        Q_rec(t,:) = agent.Q_t;                                            
        
        % if BB-H-C is evaluated
        if isequal(a_class, @th_agentobj_BB_H_C) 
            Q_rec_BB_R_1(t,:) = agent.Q_t_BB_R_1;
            Q_rec_BB_I(t,:)   = agent.Q_t_BB_I;
        % if BB-H-E is evaluated
        elseif isequal(a_class, @th_agentobj_BB_H_E) 
            Q_rec_BB_R_2(t,:) = agent.Q_t_BB_R_2;
            Q_rec_BB_I(t,:)   = agent.Q_t_BB_I;    
        end
        
        % record agent internal decision variables: agent action probabilities
        p_rec(t,:) = agent.p_t;
        
        % if BB-H-C is evaluated
        if isequal(a_class, @th_agentobj_BB_H_C) 
            p_rec_BB_R_1(t,:) = agent.p_t_BB_R_1;
            p_rec_BB_I(t,:)   = agent.p_t_BB_I;
        % if BB-H-E is evaluated
        elseif isequal(a_class, @th_agentobj_BB_H_E) 
            p_rec_BB_R_2(t,:) = agent.p_t_BB_R_2;
            p_rec_BB_I(t,:)   = agent.p_t_BB_I;
        end
        
        % record agent action
        a_rec(t) = agent.a_t; 
        
        % update trial-dependent variables
        % -----------------------------------------------------------------
        task = state_update(task, agent.a_t);                               % action dependent task state update 
        t    = t + 1;                                                       % update trial index for task and agent
              
    end
    
    % ******** evaluate the outcome of the last decision on attempt *******
    
    % task state evaluation after last decision
    % ---------------------------------------------------------------------
    ft = 0;                                                                 % initialize trial target found index
    if task.s_t(1) == task.s_t(3)
         task.fnd(1) = 1;                                                   % target 1 found
         task.tkn    = task.s_t(3);                                         % set known target location
         ft          = 1;                                                   % turn on trial target found index
         td{att}     = [td{att}, task.s_t(3)];                              % set visited target location 
    elseif task.s_t(1) == task.s_t(4)
         task.fnd(2) = 1;                                                   % target 2 found
         task.tkn    = task.s_t(4);                                         % set known target location
         ft          = 1;                                                   % turn on trial target found index
         td{att}     = [td{att}, task.s_t(4)];                              % set visited target location 
    end

    % complete final task state 
    if task.fnd(1) && ~task.fnd(2)
        task.s_t(2) = 1;                                                    % first treasure found   
    elseif task.fnd(2) && ~task.fnd(1)
        task.s_t(2) = 2;                                                    % second treasure found
    elseif task.fnd(1) && task.fnd(2)
        task.s_t(2) = 3;                                                    % both treasures found
    elseif ~task.fnd(1) && ~task.fnd(2)
        task.s_t(2) = 0;                                                    % no treasure found
    end
    
    % record state
    s_rec(t,:) = task.s_t;
    
    % agent belief state update after last decision
    % ---------------------------------------------------------------------
    agent.t   = t;                                                          % trial counter
    agent.s_t = task.s_t(1:2);                                              % observable part of the state
    agent.ft  = ft;                                                         % trial target found index
    agent.td  = td;                                                         % visited target locations 
    agent.att = task.att;                                                   % attempt on task
    agent.ap  = 2;                                                          % post trial attempt phase
    agent.tkn = task.tkn;                                                   % known target location
    
    % update belief state
    agent = belief_update(agent); 

    % record post-attempt belief state
    if ~isempty(agent.B)
        B_rec(t+1,:) = agent.B(:,end)'; 
    end
    
    % ********************* concatenate and save data *********************
    
    % unobservable data
    % ---------------------------------------------------------------------
    
    % if agent BB-H-C is evaluated
    if isequal(a_class, @th_agentobj_BB_H_C)
        x{att} = [s_rec(:,3:4) p_rec, p_rec_BB_R_1, p_rec_BB_I, Q_rec_BB_R_1, Q_rec_BB_I]; 
    % if agent BB-H-E is evaluated
    elseif isequal(a_class, @th_agentobj_BB_H_E)
        x{att} = [s_rec(:,3:4) p_rec, p_rec_BB_R_2, p_rec_BB_I, Q_rec_BB_R_2, Q_rec_BB_I]; 
    % for all other agents
    else
        x{att} = [s_rec(:,3:4) p_rec];                                      
    end
    
    % observable data
    % ---------------------------------------------------------------------
    y{att} = [s_rec(:,1:2) o_rec a_rec ones(T+1,1)*T_a(att)];              
    
    % belief state
    % ---------------------------------------------------------------------
    B{att} = B_rec;                                                         
    
    % ************** if the task is solved, end attempt loop **************
    
    if sum(task.fnd) == 2
        task.sol = 1;
        break
    end
    
end 

% ************************** set output structure *************************

a_sim.y   = y;                                                              % observable data 
a_sim.x   = x;                                                              % unobservable data
a_sim.B   = B;                                                              % belief state data
a_sim.sol = task.sol;                                                       % task solved flag

end


