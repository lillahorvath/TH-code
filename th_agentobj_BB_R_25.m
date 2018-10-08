classdef th_agentobj_BB_R_25

    % This class definition file specifies the properties and methods of an
    % object that implements agent variant BB-R-25.
    %
    % Copyright (C) Dirk Ostwald, Lilla Horvath
    % ---------------------------------------------------------------------
    
    % properties of agent variant BB-R-25
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
    
    % methods of agent variant BB-R-25
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function agentobj = th_agentobj_BB_R_25(agent_init)

        % This function creates an agent object of class 
        % th_agentobj_BB_R_25 based on the agent initialization input 
        % structure
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
        
        % This function updates the belief state of the agent object based 
        % on its current state, including a novel observation (collection 
        % of light and dark grey bars).
        %
        % Inputs
        %       agentobj: agent object with prior belief state
        %
        % Outputs
        %       agentobj: agent object with posterior belief state
        %
        % Copyright (C) Lilla Horvath, Dirk Ostwald
        % -----------------------------------------------------------------
        
        % -----------------------------------------------------------------
        %                      initial attempt phase 
        % -----------------------------------------------------------------
        
        if agentobj.ap == 0
            
            % -------------- evaluate the prior belief state --------------
            
            % first attempt on task
            % -------------------------------------------------------------
            if agentobj.att == 1

                % initialize
                B_pri = zeros(agentobj.nn,1);                              
                
                % set uniform values
                for i = 2:agentobj.nn
                    B_pri(i) = 1;                                                                            
                end
                
                % normalize
                B_pri = B_pri/sum(B_pri(:));
                
            % second or third attempt on task
            % -------------------------------------------------------------
            else
                
                % if no known target from previous attempt(s)
                % ---------------------------------------------------------
                if isnan(agentobj.tkn)
                    
                    % get the posterior belief state of the previous attempt
                    B_pri = agentobj.B(:,end); 
                    
                % if target known from previous attempt(s) 
                % ---------------------------------------------------------
                else
                    
                    % get the posterior belief state of the previous attempt
                    B_pri = agentobj.B(:,end);                
                    
                    % set known target location to maximal belief
                    B_pri(agentobj.tkn) = 1;                                
                    
                    % normalize
                    B_pri = B_pri/sum(B_pri);                 
                end
                
            end
            
            % ------------ concatenate belief state evolution -------------
            agentobj.B = [agentobj.B B_pri];
                 
        % -----------------------------------------------------------------
        %                       trial attempt phase 
        % -----------------------------------------------------------------
        
        elseif agentobj.ap == 1

            % no target found in attempt up to trial
            % -------------------------------------------------------------
            if agentobj.s_t(2) == 0 
                
                %  ----------- evaluate the prior belief state ------------
                % (1) project 1D prior belief state into 2D belief state 
                % (2) set the diagonal to zero (to reflect the agent's belief
                % that the targets are located at different nodes)
                
                % if no known target from previous attempt(s)
                % ---------------------------------------------------------
                if isnan(agentobj.tkn)
                    
                    % initialize
                    B_pri    = agentobj.B(:,end);                           
                    B_pri_2D = NaN(length(B_pri), length(B_pri));           

                    % cycle over target location pairs
                    for i = 1:length(B_pri)
                        for j = 1:length(B_pri)

                            % evaluate the joint beliefs
                            B_pri_2D(i,j) = B_pri(i)*B_pri(j);              

                            % set the diagonal to 0
                            if i == j 
                                B_pri_2D(i,j) = 0;                          
                            end

                        end
                    end
                    
                    % normalize
                    B_pri_2D = B_pri_2D/sum(B_pri_2D(:));
                    
                % if target known from previous attempt(s)
                % ---------------------------------------------------------
                else
                    
                    % initialize
                    B_pri    = agentobj.B(:,end);                           
                    B_pri_2D = zeros(length(B_pri), length(B_pri));        
                    
                    % directly assign the prior belief state
                    B_pri_2D(:,agentobj.tkn) = B_pri;                       % to column (index corresponds to the known target location)
                    B_pri_2D(agentobj.tkn,:) = B_pri';                      % to row (index corresponds to the known target location)
                    
                    % set the diagonal to 0
                    B_pri_2D(1:length(B_pri)+1:end) = 0;                    

                    % normalize
                    B_pri_2D = B_pri_2D/sum(B_pri_2D(:));                  
                   
                end
                
                % --------------- position dependent update ---------------
                
                % set the belief at the current position to zero
                B_pri_2D(agentobj.s_t(1),:) = 0;
                B_pri_2D(:,agentobj.s_t(1)) = 0;
                
                % normalize
                B_pri_2D = B_pri_2D/sum(B_pri_2D(:));
                
                % ------------- observation dependent update --------------

                % initialize
                p_o_giv_s = squeeze(agentobj.O{1}(agentobj.s_t(1),:,:,:));  % p(o_i = 1|s) (i = 1,2,3,4), NaN for not available observations
                p_s_giv_o = zeros(size(B_pri_2D));                          % p(s|[o_1,o_2,o_3,o_4]) initialization

                % cycle over target location pairs
                for l_t1 = 1:agentobj.nn
                    for l_t2 = 1:agentobj.nn

                        % determine emission probabilities based on observation
                        p_oi_giv_s = NaN(1,size(p_o_giv_s,3));

                        % cycle over observations and recover probabilities p(o_i|s)
                        for i = 1:size(p_o_giv_s,3)
                           switch agentobj.o_t(i)
                                case 1
                                    p_oi_giv_s(i) = p_o_giv_s(l_t1,l_t2,i);     % light bar observation
                                case 0                                     
                                    p_oi_giv_s(i) = 1 - p_o_giv_s(l_t1,l_t2,i); % dark bar observation
                            end
                        end

                        % change possible NaN values to 1 for use of prod
                        p_oi_giv_s(isnan(p_oi_giv_s)) = 1;

                        % evaluate Bayes theorem numerator
                        p_s_giv_o(l_t1,l_t2) = prod(p_oi_giv_s)*B_pri_2D(l_t1,l_t2);
                    end
                end

                % normalize
                p_s_giv_o = p_s_giv_o/sum(p_s_giv_o(:));

                % marginalize
                B_pos = sum(p_s_giv_o,2);

                % ---------- concatenate belief state evolution -----------
                agentobj.B = [agentobj.B B_pos];
                    
            % target found in attempt
            % -------------------------------------------------------------
            else
                
                % if target found on current trial
                % ---------------------------------------------------------
                if agentobj.ft && numel(agentobj.td{agentobj.att}) == 1 
                    
                    %  --------- evaluate the prior belief state ----------
                    % (1) project 1D prior belief state into 2D belief state 
                    % (2) set the diagonal to zero (to reflect the agent's 
                    % belief that the targets are located at different nodes)

                    % if no target known from previous attempt(s)
                    % -----------------------------------------------------
                    if ~ismember(0,isnan(agentobj.T_a))

                        % initialize
                        B_pri    = agentobj.B(:,end); 
                        B_pri_2D = NaN(length(B_pri), length(B_pri));

                        % cycle over target location pairs
                        for i = 1:length(B_pri)
                            for j = 1:length(B_pri)

                                % evaluate the joint beliefs
                                B_pri_2D(i,j) = B_pri(i)*B_pri(j);

                                % set the diagonal to 0
                                if i == j
                                    B_pri_2D(i,j) = 0;  
                                end

                            end
                        end

                        % normalize
                        B_pri_2D = B_pri_2D/sum(B_pri_2D(:));
                        
                    % if target known from a previous attempt(s)
                    % -----------------------------------------------------
                    else
                        
                        % initialize
                        B_pri    = agentobj.B(:,end); 
                        B_pri_2D = zeros(length(B_pri), length(B_pri));

                        % directly assign the prior belief state
                        B_pri_2D(:,agentobj.tkn) = B_pri;                   % to column (index corresponds to the known target location)
                        B_pri_2D(agentobj.tkn,:) = B_pri';                  % to row (index corresponds to the known target location)

                        % set the diagonal to 0
                        B_pri_2D(1:length(B_pri)+1:end) = 0;

                        % normalize
                        B_pri_2D = B_pri_2D/sum(B_pri_2D(:)); 

                    end
                    
                    % ------------- position dependent update -------------
                    
                    % marginalize the belief state over the discovered target node
                    B_pri = [];
                    B_pri = (B_pri_2D(agentobj.s_t(1),:)/sum(B_pri_2D(agentobj.s_t(1),:)))';
                    
                % if target found on a previous trial
                % ---------------------------------------------------------
                else 

                    % ---------- evaluate the prior belief state ----------
                    B_pri = agentobj.B(:,end); 

                    % ------------- position dependent update -------------

                    % set the belief at the current position to zero
                    B_pri(agentobj.s_t(1)) = 0;
                    
                    % normalize
                    B_pri = B_pri/sum(B_pri(:));
                    
                end    
                
                % -------------- observation dependent update -------------
                
                % initialize
                p_o_giv_s = squeeze(agentobj.O{2}(agentobj.s_t(1),:,:));    % p(o_i = 1|s) (i = 1,2,3,4), NaN for not available observations
                p_s_giv_o = zeros(size(B_pri));                             % p(s|[o_1,o_2,o_3,o_4]) initialization
                
                % cycle over target locations
                for l_t = 1:agentobj.nn

                    % determine emission probabilities based on observation
                    p_oi_giv_s = NaN(1,size(p_o_giv_s,2));

                    % cycle over observations and recover probabilities p(o_i|s)
                    for i = 1:size(p_o_giv_s,2)
                       switch agentobj.o_t(i)
                            case 1
                                p_oi_giv_s(i) = p_o_giv_s(l_t,i);           % light bar observation
                            case 0
                                p_oi_giv_s(i) = 1 - p_o_giv_s(l_t,i);       % dark bar observation
                        end
                    end

                    % change possible NaN values to 1 for use of prod
                    p_oi_giv_s(isnan(p_oi_giv_s)) = 1;

                    % evaluate Bayes theorem numerator
                    p_s_giv_o(l_t) = prod(p_oi_giv_s)*B_pri(l_t);

                end

                % normalize
                B_pos = p_s_giv_o/sum(p_s_giv_o(:));

                % ----------- concatenate belief state evolution ----------
                agentobj.B = [agentobj.B B_pos];

            end
            
        % -----------------------------------------------------------------
        %                     post-trial attempt phase 	
        % -----------------------------------------------------------------    
              
        elseif agentobj.ap == 2
            
            % no target found in attempt
            % -------------------------------------------------------------
            if agentobj.s_t(2) == 0 
                
               %  ------------ evaluate the prior belief state ------------ 
               % (1) project 1D prior belief state into 2D belief state 
               % (2) set the diagonal to zero (to reflect the agent's belief
               % that the targets are located at different nodes)
                
               % if no target known from previous attempt(s)
               % ----------------------------------------------------------
                if isnan(agentobj.tkn)

                    % initialize
                    B_pri    = agentobj.B(:,end); 
                    B_pri_2D = NaN(length(B_pri), length(B_pri));

                    % cycle over target location pairs
                    for i = 1:length(B_pri)
                        for j = 1:length(B_pri)

                            % evaluate the joint beliefs
                            B_pri_2D(i,j) = B_pri(i)*B_pri(j);

                            % set the diagonal to 0
                            if i == j
                                B_pri_2D(i,j) = 0;  
                            end

                        end
                    end

                    % normalize
                    B_pri_2D = B_pri_2D/sum(B_pri_2D(:));
                    
                % if target known from previous attempt(s)
                % ---------------------------------------------------------
                else
                    
                    % initialize
                    B_pri    = agentobj.B(:,end);
                    B_pri_2D = zeros(length(B_pri), length(B_pri));
                    
                    % directly assign the prior belief state
                    B_pri_2D(:,agentobj.tkn) = B_pri;                       % to column (index corresponds to the known target location)
                    B_pri_2D(agentobj.tkn,:) = B_pri';                      % to row (index corresponds to the known target location)
                    
                    % set the diagonal to 0
                    B_pri_2D(1:length(B_pri)+1:end) = 0;
                    
                    % normalize
                    B_pri_2D = B_pri_2D/sum(B_pri_2D(:));
                   
                end
                                
                % --------------- position dependent update ---------------
                
                % set the belief at the current position to zero
                B_pri_2D(agentobj.s_t(1),:) = 0;
                B_pri_2D(:,agentobj.s_t(1)) = 0;
                
                % normalize
                B_pri_2D = B_pri_2D/sum(B_pri_2D(:));
                
                % marginalize
                B_pos = sum(B_pri_2D,2);

                % ----------- concatenate belief state evolution ----------
                agentobj.B = [agentobj.B B_pos];
                
            % one target found
            % -------------------------------------------------------------
            elseif agentobj.s_t(2) == 1 || agentobj.s_t(2) == 2

                % if target found on current trial
                % ---------------------------------------------------------
                if agentobj.ft && numel(agentobj.td{agentobj.att}) == 1 
                    
                    %  --------- evaluate the prior belief state ---------- 
                    % (1) project 1D prior belief state into 2D belief state 
                    % (2) set the diagonal to zero (to reflect the agent's 
                    % belief that the targets are located at different nodes)
                    
                    % if no target known from previous attempt(s)
                    % -----------------------------------------------------
                    if ~ismember(0,isnan(agentobj.T_a)) 

                        % initialize
                        B_pri    = agentobj.B(:,end); 
                        B_pri_2D = NaN(length(B_pri), length(B_pri));

                        % cycle over target location pairs
                        for i = 1:length(B_pri)
                            for j = 1:length(B_pri)

                                % evaluate the joint beliefs
                                B_pri_2D(i,j) = B_pri(i)*B_pri(j);

                                % set the diagonal to 0
                                if i == j
                                    B_pri_2D(i,j) = 0;  
                                end

                            end
                        end

                        % normalize
                        B_pri_2D = B_pri_2D/sum(B_pri_2D(:));
                    
                    %  if target known from previous attempt(s)
                    % -----------------------------------------------------
                    else
                        
                        % initialize
                        B_pri    = agentobj.B(:,end);
                        B_pri_2D = zeros(length(B_pri), length(B_pri));

                        % directly assign the prior belief state
                        B_pri_2D(:,agentobj.tkn) = B_pri;                   % to column (index corresponds to the known target location)
                        B_pri_2D(agentobj.tkn,:) = B_pri';                  % to row (index corresponds to the known target location)
                        
                        % set the diagonal to 0
                        B_pri_2D(1:length(B_pri)+1:end) = 0;

                        % normalize
                        B_pri_2D = B_pri_2D/sum(B_pri_2D(:)); 

                    end
                    
                    % ------------- position dependent update -------------
                    
                    % marginalize the belief state over the discovered target node
                    B_pri = [];
                    B_pri = (B_pri_2D(agentobj.s_t(1),:)/sum(B_pri_2D(agentobj.s_t(1),:)))';
                    
                    % --------- concatenate belief state evolution --------
                    B_pos      = B_pri;
                    agentobj.B = [agentobj.B B_pos];
                    
                % if target found on a previous trial
                % ---------------------------------------------------------
                else 
                    
                    %  --------- evaluate the prior belief state ---------- 
                    B_pri = agentobj.B(:,end); 

                    % ------------- position dependent update -------------
                    
                    % set the belief at the current position to zero
                    B_pri(agentobj.s_t(1)) = 0;
                    
                    % normalize
                    B_pri = B_pri/sum(B_pri(:));
                    
                    % --------- concatenate belief state evolution --------
                    B_pos      = B_pri;
                    agentobj.B = [agentobj.B B_pos];

                end
                
            % if task solved 
            % -------------------------------------------------------------
            elseif agentobj.s_t(2) == 3
                
                % ------------ evaluate the prior belief state ------------
                B_pri = agentobj.B(:,end); 
                
                % --------------- position dependent update ---------------
                
                % set the belief at the current position to one and to zero 
                % for all other nodes
                B_pos                  = zeros(length(B_pri),1);
                B_pos(agentobj.s_t(1)) = 1;

                % ---------- concatenate belief state evolution -----------
                agentobj.B = [agentobj.B B_pos];
                
            end
            
        end
        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function agentobj = decide(agentobj)
 
        % This function implements the agent's decision strategy based its 
        % current belief state. It uses a one-step look ahead procedure 
        % based on the complete belief state distribution to evaluate the
        % action valences
        %
        % Inputs
        %       agentobj: pre-decision agent object
        %
        % Outputs
        %       agentobj: decided agent object
        %
        % Copyright (C) Dirk Ostwald, Lilla Horvath
        % -----------------------------------------------------------------
        
        % initialization
        % -----------------------------------------------------------------
        A_s          = find(~isnan(agentobj.o_t));                          % indices of available actions
        V_a_s        = NaN(1,length(agentobj.A));                           % state-action valence
        agentobj.Q_t = NaN(1,length(agentobj.A));                           % action valence
        agentobj.p_t = NaN(1,length(agentobj.A));                           % action probability
        
        % action valence evaluation: compute the belief state weighted l1
        % distances 
        % -----------------------------------------------------------------
        % cycle over available actions
        for i = 1:length(A_s)

            % one-step look ahead
            s_t_t = agentobj.s_t(1) + agentobj.A(A_s(i));

            % target distances
            V_a_s(A_s(i)) = agentobj.L1(s_t_t,:)*round(agentobj.B(:,end), 10);

        end

        % set action valences
        agentobj.Q_t = -V_a_s;
  
        % action probability evaluation
        % -----------------------------------------------------------------
        agentobj.p_t(A_s) = th_softmax(agentobj.Q_t(A_s),agentobj.be);   

        % action sampling from available action set
        % -----------------------------------------------------------------
        agentobj.a_t = agentobj.A(A_s(logical(mnrnd(1,agentobj.p_t(A_s)))));
      
        end
    end
end
        
        