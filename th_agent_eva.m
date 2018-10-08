function a_eva = th_agent_eva(a_eva)

% This function evaluates the choice data of the attempts of a single task
% of the treasure hunt task with respect to its likelihood under an agent 
% model M and records internal agent variables for parametric GLM-fMRI
% analyses.
%
% Inputs
%       a_eva:      structure with fields
%        .pomdp:    POMDP formulation
%        .y:        observable data (observed state-observation-action evolution)
%        .a_class:  agent class
%        .par:      parameter values 
%
% Outputs
%       a_eva:      input structure with additional fields
%        .ll:       data set log likelihood
%        .nva:      number of valid actions
%        .B_evo:    belief state evolution
%        .Q_evo:    chosen action valence evolution
%        .NAC_evo:  normalized choice conflict evolution
%
% Copyright (C) Dirk Ostwald, Lilla Horvath
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%                              initialization 
% -------------------------------------------------------------------------

% unpack input structure
y       = a_eva.y;                                                          
a_class = a_eva.a_class;                                                    
A       = a_eva.pomdp.A;                                                    % action set
A_s     = a_eva.pomdp.A_s;                                                  % position dependent available action set
natt    = length(find(~cellfun(@isempty,y)));                               % number of non-empty task attempts
lpe     = NaN(natt,1);                                                      % log probability array

% initialize task-attempt arrays and counters
B_evo   = cell(natt,1);                                                     % belief state evolution
Q_evo   = cell(natt,1);                                                     % chosen action valence evolution
NAC_evo = cell(natt,1);                                                     % normalized action entropy evolution
nva     = 0;                                                                % number of valid actions counter
td      = cell(natt,1);                                                     % cell array with target location visits per attempts

% ************************ initialize agent object ************************

% specify agent initialization structure
a_init.nn  = a_eva.pomdp.nn;                                                % number of gridworld nodes
a_init.L1  = a_eva.pomdp.L1;                                                % L1 distances
a_init.A   = a_eva.pomdp.A;                                                 % set of actions
a_init.O   = a_eva.pomdp.O;                                                 % observation likelihood function
a_init.att = 0;                                                             % initialize attempt counter
a_init.B   = [];                                                            % initialize belief state
a_init.tkn = NaN;                                                           % initialize known target location
a_init.t   = NaN;                                                           % initialize trial counter

if isequal(a_class, @th_agentobj_BF_0)
    a_init.be = a_eva.par;                                                  % BF-0: no parameter 
else
    if isequal(a_class, @th_agentobj_BB_H_C)
        a_init.la = a_eva.par(2);                                           % BB-H-C: lambda
    elseif isequal(a_class, @th_agentobj_BB_H_E)     
        a_init.os = a_eva.par(2);                                           % BB-H-E: offset
        a_init.sl = a_eva.par(3);                                           % BB-H-E: slope
    end    
    a_init.be = a_eva.par(1);                                               % all agents except BF-0: beta
end

% initialize agent
agent = feval(a_class, a_init);                                             % agent object creation

% -------------------------------------------------------------------------
%                              evaluate agent 
% -------------------------------------------------------------------------

% cycle over attempts
for att = 1:natt
    
    % ******************* initialize attempt variables ********************
    
    yatt          = y{att}(1:end,:);                                        % select all standard trials    
    T_a(att)      = yatt(1,8);                                              % known target location from previous attempt(s)
    lp            = NaN(size(yatt,1)-1,1);                                  % chosen action log probability array 
    B_evo_a       = cell(size(yatt,1)+1,1);                                 % belief state attempt evolution
    act_confl_a_n = NaN(size(yatt,1)-1,1);                                  % normalized action entropy attempt evolution
    Q_evo_a       = NaN(size(yatt,1)-1,1);                                  % chosen action valence attempt evolution
    
    f_v  = find(logical(yatt(1:end,2)~=0),1);                               % trial number of first target location visit in attempt
    if ~isempty(f_v)
        a_v = find(yatt(:,1) == yatt(f_v,1));                               % list of trial numbers where target is visited
    else
        a_v = 0;
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
    if ~isempty(agent.B)
        B_evo_a{1} = agent.B(:,end);                                        % attempt-specific prior belief state
    end
    
    % cycle over trials 
    for t = 1:size(yatt,1)-1
        
        % ************ evaluate choice log likelihood on trial ************
        
        % set key trial data
        % -----------------------------------------------------------------
        s_t = yatt(t,1:2);                                                  % observable part of the state: position and attempt target found flag
        o_t = yatt(t,3:6);                                                  % observation bars
        a_t = yatt(t,7);                                                    % action
        tkn = yatt(t,8);                                                    % known target location
        
        ft = 0;                                                             % initialize trial target found index
        if ismember(t,a_v)
            ft      = 1;                                                    % turn on trial target found index
            td{att} = [td{att}, s_t(1)];                                    % set visited target location 
        end
        
        % evaluate agent dependent action probabilities
        % -----------------------------------------------------------------
        agent.t   = t;                                                      % trial counter
        agent.s_t = s_t;                                                    % observable part of the state: position and attempt target found flag
        agent.ft  = ft;                                                     % trial target found index
        agent.td  = td;                                                     % visited target locations
        agent.o_t = o_t;                                                    % observation bars
        agent.att = att;                                                    % attempt on task
        agent.ap  = 1;                                                      % trial attempt phase
        agent.tkn = tkn;                                                    % known target location 
        
        % evaluate belief state
        agent = belief_update(agent); 
        
        % record updated belief state, if existent
        if ~isempty(agent.B)
            B_evo_a{t+1} = agent.B(:,end); 
        end
        
        % evaluate action probabilities
        agent = decide(agent); 
           
        % for valid responses
        if ~isnan(a_t) && ismember(a_t,A(logical(A_s(s_t(1),:))))
            
            % evaluate log probability of the chosen action
            lp(t) = log(agent.p_t(A == a_t));
            
            % get normalized action conflict: measured as the normalized entropy of the softmax distibuted action RV
            act_confl_a_n(t) = -nansum(agent.p_t.*log(agent.p_t)./log(sum(~isnan(agent.p_t))));
        
            % record chosen action valence
            Q_evo_a(t) = agent.Q_t(A == a_t);
            
            % update valid actions counter
            nva = nva + 1;

        % if no, or an invalid response was recorded, set ln p(a_t) = 0
        else
            lp(t) = 0;
        end

    end
    
    % ******** evaluate the outcome of the last decision on attempt *******
    
    % set key last trial data
    % ---------------------------------------------------------------------
    t   = t+1;                                                              % trial counter
    s_t = yatt(t,1:2);                                                      % observable part of the state: position and attempt target found flag
    tkn = yatt(t,8);                                                        % known target location

    ft = 0;                                                                 % initialize trial target found index
    if ismember(t,a_v)
        ft      = 1;                                                        % turn on trial target found index
        td{att} = [td{att}, s_t(1)];                                        % set visited target location
    end
    
    % belief state update after last decision
    % ---------------------------------------------------------------------
    agent.t   = t;                                                          % trial counter
    agent.s_t = s_t;                                                        % observable part of the state: position and attempt target found flag
    agent.ap  = 2;                                                          % post trial attempt phase
    agent.tkn = tkn;                                                        % known target node
    agent.ft  = ft;                                                         % trial target found index
    agent.td  = td;                                                         % visited target locations
    agent.att = att;                                                        % attempt on task
    
    % update belief state
    agent = belief_update(agent);                                           

    % record post-attempt belief state
    if ~isempty(agent.B)
        B_evo_a{t+1} = agent.B(:,end); 
    end
    
    % ************************* save attempt data *************************
    
    B_evo{att}   = B_evo_a;
    Q_evo{att}   = Q_evo_a;
    NAC_evo{att} = act_confl_a_n;
      
    % evaluate choice data log likelihood for the current attempt
    lpe(att) = sum(lp);

end

% ************************** set output structure *************************

a_eva.ll      = sum(lpe);                                                   % task choice data log likelihood
a_eva.nva     = nva;                                                        % number of valid responses
a_eva.B_evo   = B_evo;                                                      % belief state evolution
a_eva.Q_evo   = Q_evo;                                                      % chasen action valence evolution
a_eva.NAC_evo = NAC_evo;                                                    % normalized choice conflict evolution

end

