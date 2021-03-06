classdef th_agentobj_BF_0

    % This class definition file specifies the properties and methods of an
    % object that implements agent variant BF-0.
    %
    % Copyright (C) Dirk Ostwald, Lilla Horvath
    % ---------------------------------------------------------------------
    
    % properties of agent variant BF-0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        
        nn;                 % number of grid-world nodes
        L1;                 % L1 distances
        A;                  % available action set
        O;                  % observation likelihood function
        be;                 % softmax parameter beta
        att;                % task attempt counter
        ap;                 % attempt phase marker
        t;                  % trial counter
        ft;                 % target found on trial index
        B;                  % belief state evolution
        T_a;                % array of known target locations from previous attempt(s) 
        td;                 % cell array with target location visits per attempts
        tkn;                % known target location (latest discovery)
        s_t;                % observable state (= position and attempt target found flag)
        o_t;                % observation state
        Q_t;                % action valence state
        p_t;                % action probability state
        a_t;                % action state
         
    end
    
    % methods of agent variant BF-0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function agentobj = th_agentobj_BF_0(agent_init)

        % This function creates an agent object of class th_agentobj_BF_0 
        % based on the agent initialization input structure.
        %
        %   Inputs
        %       agent_init: structure with agent initialization fields
        %
        %   Outputs
        %       agentobj:   agent object with initialized values
        %
        % Copyright (C) Dirk Ostwald, Lilla Horvath
        % -----------------------------------------------------------------      
        
        % set agent properties based on input structure
        agentobj.nn  = agent_init.nn;
        agentobj.L1  = agent_init.L1;
        agentobj.A   = agent_init.A;
        agentobj.O   = agent_init.O;
        agentobj.be  = agent_init.be;
        agentobj.att = agent_init.att;
        agentobj.B   = agent_init.B;
        agentobj.tkn = agent_init.tkn;
       
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function agentobj = belief_update(agentobj)
        
        % Agent BF-0 does not have a belief state, the function is thus
        % empty.
        %
        %   Inputs
        %       agentobj:   agent object without belief state
        %
        %   Outputs
        %       agentobj:   agent object without belief state
        %
        % Copyright (C) Dirk Ostwald, Lilla Horvath
        % -----------------------------------------------------------------
        return
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function agentobj = decide(agentobj)
 
        % This function implements the agent's decision strategy of choosing
        % actions at random.
        %
        %   Inputs
        %       agentobj:   pre-decision agent object
        %
        %   Outputs
        %       agentobj:   decided agent object
        %
        % Copyright (C) Dirk Ostwald, Lilla Horvath
        % -----------------------------------------------------------------
        
        % initialization
        % -----------------------------------------------------------------
        A_s          = find(~isnan(agentobj.o_t));                          % indicies of available actions
        agentobj.Q_t = NaN(1,length(agentobj.A));                           % action valences
        
        % action valence evaluation: assign equal valences to available actions
        % -----------------------------------------------------------------
        agentobj.Q_t(A_s) = 1/length(A_s);

        % action probability evaluation
        % -----------------------------------------------------------------
        agentobj.p_t = agentobj.Q_t;
        
        % action sampling from available action set
        % -----------------------------------------------------------------
        agentobj.a_t = agentobj.A(A_s(logical(mnrnd(1,agentobj.p_t(A_s)))));
   
        end
    end
end
        
        