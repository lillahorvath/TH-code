classdef th_agentobj_BB_H_E

    % This class definition file specifies the properties and methods of an
    % object that implements agent variant BB-H-E.
    %
    % Copyright (C) Dirk Ostwald, Lilla Horvath
    % ---------------------------------------------------------------------
    
    % properties of agent variant BB-H-E
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        
        nn;                 % number of grid-world nodes
        L1;                 % L1 distances
        A;                  % available action set
        O;                  % observation likelihood function
        be;                 % softmax parameter beta
        os;                 % exponential decay offset parameter 
        sl;                 % exponential decay slope parameter
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
        Q_t;                % action valence state (hybrid)
        Q_t_BB_R_2;         % action valence (component: BB-R-1)
        Q_t_BB_I;           % action valence (component: BB-I)
        p_t;                % action probability state (hybrid)
        p_t_BB_R_2;         % action probability (component: BB-R-1)
        p_t_BB_I;           % action probability (component: BB-I)
        a_t;                % action state
        s_t_vb;             % virtually simulated one-step lookahead observable state
        o_t_vb;             % virtually simulated one-step lookahead observation   
         
    end
    
    % methods of agent variant BB-H-E
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function agentobj = th_agentobj_BB_H_E(agent_init)

        % This function creates an agent object of class th_agentobj_BB_H_E        
        % based on the agent initialization input structure.
        %
        % Inputs
        %       agent_init: structure with agent initialization fields
        %
        % Outputs
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
        agentobj.os  = agent_init.os;
        agentobj.sl  = agent_init.sl;
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
 
        % This function implements the agent's decision strategy, which is
        % a hybrid strategy combining the actions' information-based (as 
        % computed with the BB-I scheme) and reward-based valences (as 
        % computed with the BB-R-2 scheme). This hybrid agent uses a state
        % uncertainty-dependent parameter that dynamically controls the 
        % linear combination of the two action valences.
        %
        % Inputs
        %       agentobj: pre-decision agent object
        %
        % Outputs
        %       agentobj: decided agent object
        %
        % Copyright (C) Lilla Horvath, Dirk Ostwald
        % -----------------------------------------------------------------        
        
        % -----------------------------------------------------------------
        %                           initialization 
        % -----------------------------------------------------------------
        A_s                 = find(~isnan(agentobj.o_t));                   % indices of available actions
        V_a_s_BB_I          = NaN(1,length(agentobj.A));                    % BB-I state-action valence
        n_V_a_s_BB_I        = NaN(1,length(agentobj.A));                    % BB-I normalized state-action valence
        V_a_s_BB_R_2        = NaN(1,length(agentobj.A));                    % BB-R-2 state-action valence
        n_V_a_s_BB_R_2      = NaN(1,length(agentobj.A));                    % BB-R-2 normalized state-action valence
        agentobj.Q_t        = NaN(1,length(agentobj.A));                    % BB-H-C action valence 
        agentobj.Q_t_BB_I   = NaN(1,length(agentobj.A));                    % BB-I action valence (normalized)
        agentobj.Q_t_BB_R_2 = NaN(1,length(agentobj.A));                    % BB-R-2 action valence (normalized)
        agentobj.p_t        = NaN(1,length(agentobj.A));                    % BB-H-C action probability
        agentobj.p_t_BB_I   = NaN(1,length(agentobj.A));                    % BB-I action probability 
        agentobj.p_t_BB_R_2 = NaN(1,length(agentobj.A));                    % BB-R-2 action probability
        curr_bf             = round(agentobj.B(:,end),10);                  % get the agent's current belief
        
        % -----------------------------------------------------------------
        %                     action valence evaluation 
        % -----------------------------------------------------------------
        
        % ************* informative valence: the BB-I scheme **************
        
        % scale the belief at a certain node depending of the attempt
        % specific target discovery history
        if agentobj.s_t(2) == 0
            mbn = 2;
        else
            mbn = 1;
        end
        
        % cycle over available actions
        for a = 1:length(A_s)
                       
            s_t_t     = agentobj.s_t(1) + agentobj.A(A_s(a));               % one-step look ahead position            
            VB_pos_pd = NaN(agentobj.nn,2);                                 % virtual position dependent belief state 

            % assume no target discovery in the next trial
            % -------------------------------------------------------------
            
            % ----------- position dependent belief state update ----------
            
            % set virtually simulated one-step lookahead observable state
            agentobj.s_t_vb = [s_t_t, 0]; 
            
            % get virtual position dependent belief state
            [VB_pos_pd(:,1), VB_pos_pd_2D_nd] = pd_virtual_belief_update(agentobj);
            vbs_pd_nd                         = VB_pos_pd(:,1);
            vbs_pd_nd                         = round(vbs_pd_nd, 10);       % round to avoid an infinite value

            % KL divergence: virtual position dependent belief state and
            % current trial posterior
            KL_pd_nd = 0;                                                   % initialize  
            for pos = 1:agentobj.nn                   
                % exclude zero probability grid nodes
                if curr_bf(pos) ~= 0 && vbs_pd_nd(pos) ~= 0
                    KL_pd_nd = KL_pd_nd + vbs_pd_nd(pos)*log(vbs_pd_nd(pos)/curr_bf(pos));
                end                           
            end
            
            % -------------- get set of possible observations -------------
            
            % initialize
            obs_val   = [1,0];                                              % possible observation bar values             
            max_obs   = numel(agentobj.o_t);                                % max number of observation bars
            max_obs_m = ones(numel(obs_val)^max_obs,max_obs);               % observation matrix with all possible bar values
            
            % max set of observation bar combinations
            c = 1;                                                          % counter
            for i = 1:numel(obs_val)
                for j = 1:numel(obs_val)
                    for k = 1:numel(obs_val)
                        for l = 1:numel(obs_val)
                           max_obs_m(c,:) = [obs_val(l) obs_val(k) obs_val(j) obs_val(i)];
                           c = c+1;
                        end
                    end
                end
            end
            
            % available observations at the virtual next position
            w_o_vb = [];
            w_o_vb = find(~isnan(squeeze(agentobj.O{1}(s_t_t,1,1,:))));
            
            % number of available observations at the virt next pos 
            num_o_vb = [];                                                     
            num_o_vb = numel(w_o_vb); 
            
            % set of observation bar combinations at the virt next pos
            obs_m_vb = [];  
            obs_m_vb = max_obs_m(1:(numel(obs_val)^num_o_vb),1:num_o_vb); 
            
            % setsize
            num_obs_config = []; 
            num_obs_config = numel(obs_val)^num_o_vb; 
            
            % --------- observation dependent belief state update ---------
            
            % initialize
            VB_pos_od_nd = NaN(agentobj.nn, num_obs_config);                % virtual observation dependent belief state
            p_o_vb_nd    = NaN(num_obs_config,1);                           % observation probability under all possible grid node position values
            KL_od_nd_all = NaN(num_obs_config,1);                           % KL divergences

            for i = 1:num_obs_config                    

                % set virtually simulated one-step lookahead observation
                obs_vb          = NaN(1,4);                                 % initialize observation bar combination
                obs_vb(w_o_vb)  = obs_m_vb(i,:);   
                agentobj.o_t_vb = obs_vb;                                   
                
                % get virtual observation dependent belief state
                [VB_pos_od_nd(:,i), p_o_vb_nd(i)] = od_virtual_belief_update(agentobj, VB_pos_pd(:,1), VB_pos_pd_2D_nd);
                vbs_od_nd                         = [];
                vbs_od_nd                         = VB_pos_od_nd(:,i);
                vbs_od_nd                         = round(vbs_od_nd, 10);

                % KL divergence: virtual observation dependent belief state
                % and virtual position dependent belief state 
                KL_od_nd = 0;
                for pos = 1:agentobj.nn                   
                    % exclude zero probability grid nodes
                    if vbs_pd_nd(pos) ~= 0 && vbs_od_nd(pos) ~= 0
                        KL_od_nd = KL_od_nd + vbs_od_nd(pos)*log(vbs_od_nd(pos)/vbs_pd_nd(pos));
                    end                           
                end
                KL_od_nd_all(i) = round(p_o_vb_nd(i), 10) * KL_od_nd;       % virtual observation bar data weighted KL divergence
            end

            KL_od_nd_sum = sum(KL_od_nd_all);                               % sum of KL divergences (obs dep vb || pd dep vb)            
            KL_nd = (1-(mbn*curr_bf(s_t_t))) * (KL_pd_nd + KL_od_nd_sum);   % weighted sum of KL divergences: no target discovery assumed 

            % assume target discovery in the next trial
            % -------------------------------------------------------------
            
            % ----------- position dependent belief state update ----------
            
            % set virtually simulated one-step lookahead observable state
            agentobj.s_t_vb = [s_t_t, 1]; 
            
            % get virtual position dependent belief state
            [VB_pos_pd(:,2), VB_pos_pd_2D_d] = pd_virtual_belief_update(agentobj);
            vbs_pd_d                         = VB_pos_pd(:,2);
            vbs_pd_d                         = round(vbs_pd_d, 10);

            % KL divergence: virtual position dependent belief state and
            % current trial posterior
            KL_pd_d = 0;                                                    % initialize  
            for pos = 1:agentobj.nn                   
                % exclude zero probability grid nodes
                if curr_bf(pos) ~= 0 && vbs_pd_d(pos) ~= 0
                    KL_pd_d = KL_pd_d + vbs_pd_d(pos)*log(vbs_pd_d(pos)/curr_bf(pos));
                end                           
            end
            
            % --------- observation dependent belief state update ---------
            
            % initialize
            VB_pos_od_d = NaN(agentobj.nn, num_obs_config);                 % virtual observation dependent belief state 
            p_o_vb_d    = NaN(num_obs_config,1);                            % observation probability under all possible location (grid node) values
            KL_od_d_all = NaN(num_obs_config,1);                            % KL divergences

            for i = 1:num_obs_config                   
                
                % set virtually simulated one-step lookahead observation
                obs_vb          = NaN(1,4);                                 % initialize observation bar combination
                obs_vb(w_o_vb)  = obs_m_vb(i,:); 
                agentobj.o_t_vb = obs_vb;
                
                % get virtual observation dependent belief state
                [VB_pos_od_d(:,i), p_o_vb_d(i)] = od_virtual_belief_update(agentobj, VB_pos_pd(:,2), VB_pos_pd_2D_d);
                vbs_od_d                        = []; 
                vbs_od_d                        = VB_pos_od_d(:,i);
                vbs_od_d                        = round(vbs_od_d, 10);

                % KL divergence: virtual observation dependent belief state
                % and virtual position dependent belief state 
                KL_od_d = 0;
                for pos = 1:agentobj.nn                   
                        % exclude zero probability grid locations
                        if vbs_pd_d(pos) ~= 0 && vbs_od_d(pos) ~= 0
                            KL_od_d = KL_od_d + vbs_od_d(pos)*log(vbs_od_d(pos)/vbs_pd_d(pos));
                        end                           
                end
                KL_od_d_all(i) = round(p_o_vb_d(i), 10) * KL_od_d;          % virtual observation bar data weighted KL divergence

            end

            KL_od_d_sum = sum(KL_od_d_all);                                 % sum of KL divergences (obs dep vb || ld dep vb)            
            KL_d = (mbn*curr_bf(s_t_t)) * (KL_pd_d + KL_od_d_sum);          % weighted sum of KL divergences: no target discovery assumed
            
            V_a_s_BB_I(A_s(a)) = KL_nd + KL_d;
            
        end
        
        % normalize
        if nansum(V_a_s_BB_I) ~= 0
            n_V_a_s_BB_I  = V_a_s_BB_I/nansum(V_a_s_BB_I);
        else 
            n_V_a_s_BB_I(A_s) = 1/numel(A_s);
        end
        
        % ************* rewarding valence: the BB-R-2 scheme **************
        
        % sort and list the index values
        [~,sort_ind] = sort(curr_bf(~isnan(curr_bf)));                   
        
        % if two remaining targets in attempt: identify the two maximum-a-
        % posteriori nodes and get the one that is closer 
        if agentobj.s_t(2) == 0 
            
            [~,min_ind] = min([agentobj.L1(agentobj.s_t(1),sort_ind(end-1)), agentobj.L1(agentobj.s_t(1),sort_ind(end))]);
            if min_ind == 1
                map_n = sort_ind(end-1);
            elseif min_ind == 2
                map_n = sort_ind(end);
            end
            
        % if one remaining target in attempt: get the maximum-a-posterriori 
        % node   
        else
            map_n = sort_ind(end);
        end

        % cycle over available actions
        for i = 1:length(A_s)

            % one-step lookahead
            s_t_t = agentobj.s_t(1) + agentobj.A(A_s(i));

            % target distances
            V_a_s_BB_R_2(A_s(i)) = agentobj.L1(s_t_t,map_n);

        end
        
        % normalize
        n_V_a_s_BB_R_2 = V_a_s_BB_R_2/nansum(V_a_s_BB_R_2);

        % ********************* hybrid action valence *********************

        % compute the lambda weighting parameter
        % -----------------------------------------------------------------

        % beleif state entropy 
        bs_entropy = [];
        bs_entropy = -nansum(curr_bf.*log(curr_bf));
        
        % set the belief state entropy to zero of both treasure locations
        % are known
        if length(find(curr_bf==0.5))==2
            bs_entropy = 0;
        end
        
        % express lambda using an exponential decay function
        lambda = agentobj.os * exp(-agentobj.sl*bs_entropy);
                
        % compute hybrid action valences 
        % -----------------------------------------------------------------

        % cycle over available actions
        for a = 1:length(A_s)            
           
            % set hybrid action valences
            agentobj.Q_t(A_s(a)) = lambda*-n_V_a_s_BB_R_2(A_s(a)) + (1-lambda)*n_V_a_s_BB_I(A_s(a));
            
        end
        
        % set normalized action valences
        agentobj.Q_t_BB_I   = n_V_a_s_BB_I;
        agentobj.Q_t_BB_R_2 = -n_V_a_s_BB_R_2;

        % -----------------------------------------------------------------
        %                   action probability evaluation 	
        % -----------------------------------------------------------------
        agentobj.p_t(A_s) = th_softmax(agentobj.Q_t(A_s),agentobj.be);           % BB-H-C
        
        agentobj.p_t_BB_I(A_s) = th_softmax(n_V_a_s_BB_I(A_s),agentobj.be);      % BB-I
        agentobj.p_t_BB_R_2(A_s) = th_softmax(-n_V_a_s_BB_R_2(A_s),agentobj.be); % BB-R-2
        
        % -----------------------------------------------------------------
        %             action sampling from available action set 
        % -----------------------------------------------------------------
        agentobj.a_t = agentobj.A(A_s(logical(mnrnd(1,agentobj.p_t(A_s)))));

        end

        function [VB_pos_pd_1D, VB_pos_pd_2D] = pd_virtual_belief_update(agentobj)
            
        % This is a helper function of the agent object's decide method. It
        % evaluates the virtually simulated one-step lookahead position
        % updated belief state. 
        %
        % Inputs
        %       agentobj:     agent object without virtually simulated 
        %                     position dependent belief state 
        %
        % Outputs
        %       VB_pos_pd_1D: 1D virtually simulated position updated 
        %                     belief state
        %       VB_pos_pd_2D: 2D virtually simulated position updated 
        %                     belief state (empty if one or no remaining 
        %                     target assumed)
        %
        % Copyright (C) Lilla Horvath, Dirk Ostwald
        % -----------------------------------------------------------------
            
        % no target found in attempt up to trial
        % -----------------------------------------------------------------
        if agentobj.s_t(2) == 0

            % if no target discovery assumed
            % -------------------------------------------------------------
            if agentobj.s_t_vb(2) == 0
                
                %  ------- evaluate the virtual prior belief state -------- 
                % (1) project 1D prior belief state into 2D belief state 
                % (2) set the diagonal to zero (to reflect the agent's belief
                % that the targets are located at different nodes)

                % if no know target from previous attempt(s)
                % ---------------------------------------------------------
                if isnan(agentobj.tkn)

                    % initialize
                    B_pri    = agentobj.B(:,end); 
                    B_pri_2D = NaN(length(B_pri), length(B_pri));

                    % cycle over target location pairs
                    for i = 1:length(B_pri)
                        for j = 1:length(B_pri)

                            % evaluate joint beliefs
                            B_pri_2D(i,j) = B_pri(i)*B_pri(j);

                            % set the diagonal to 0
                            if i == j 
                                B_pri_2D(i,j) = 0;  
                            end

                        end
                    end

                    % normalize
                    if sum(B_pri_2D(:)) ~= 0
                        B_pri_2D = B_pri_2D/sum(B_pri_2D(:));
                    end

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
                    if sum(B_pri_2D(:)) ~= 0
                        B_pri_2D = B_pri_2D/sum(B_pri_2D(:));
                    end

                end

                % ----------- virtual position dependent update -----------
                
                % set the belief at the virtual position to zero
                B_pri_2D(agentobj.s_t_vb(1),:) = 0;
                B_pri_2D(:,agentobj.s_t_vb(1)) = 0;
                
                % normalize
                if sum(B_pri_2D(:)) ~= 0
                    B_pri_2D = B_pri_2D/sum(B_pri_2D(:));
                else
                    B_pri_2D = zeros(agentobj.nn, agentobj.nn);
                end
                
                % set 2D posterior
                VB_pos_pd_2D = B_pri_2D;

                % marginalize
                VB_pos_pd_1D = sum(VB_pos_pd_2D,2);

            % if target discovery assumed
            % -------------------------------------------------------------
            elseif agentobj.s_t_vb(2) == 1
                
                % -------- evaluate the virtual prior belief state -------- 
                % (1) project 1D prior belief state into 2D belief state 
                % (2) set the diagonal to zero (to reflect the agent's belief
                % that the targets are located at different nodes)

                % if no known target from previous attempt(s)
                % ---------------------------------------------------------
                if ~ismember(0,isnan(agentobj.T_a)) 

                    %  initialize
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
                    if sum(B_pri_2D(:)) ~= 0
                        B_pri_2D = B_pri_2D/sum(B_pri_2D(:));
                    end

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
                    if sum(B_pri_2D(:)) ~= 0
                        B_pri_2D = B_pri_2D/sum(B_pri_2D(:));
                    end

                end
                
                % ----------- virtual position dependent update -----------
                
                % marginalize the belief state over the virtually 
                % discovered target node
                if sum(B_pri_2D(agentobj.s_t_vb(1),:)) ~= 0
                    VB_pos_pd_1D = (B_pri_2D(agentobj.s_t_vb(1),:)/sum(B_pri_2D(agentobj.s_t_vb(1),:)))';
                else 
                    VB_pos_pd_1D = zeros(agentobj.nn,1);
                end
                
                % set 2D posterior
                VB_pos_pd_2D = [];

            end
            
        % target found in attempt
        % -----------------------------------------------------------------    
        else
            
            % if no target discovery assumed
            % ------------------------------------------------------------- 
            if agentobj.s_t_vb(2) == 0

                % -------- evaluate the virtual prior belief state --------
                B_pri = agentobj.B(:,end); 

                % ----------- virtual position dependent update ----------- 
                
                % set the belief at the virtual position to zero
                B_pri(agentobj.s_t_vb(1)) = 0;
                
                % normalize
                if sum(B_pri(:)) ~= 0
                    VB_pos_pd_1D = B_pri/sum(B_pri(:));
                else
                    VB_pos_pd_1D = zeros(agentobj.nn,1);
                end
                
                % set 2D posterior
                VB_pos_pd_2D = [];

            % if target discovery assumed
            % -------------------------------------------------------------
            elseif agentobj.s_t_vb(2) == 1 

                % ----------- virtual position dependent update -----------
                
                % set the belief at the virtual position to one and to zero 
                % for all other nodes
                VB_pos_pd_1D                     = zeros(agentobj.nn,1);
                VB_pos_pd_1D(agentobj.s_t_vb(1)) = 1;
                
                % set 2D posterior
                VB_pos_pd_2D = [];
                
            end

        end
        end
        
        function [VB_pos_od, p_o_vb] = od_virtual_belief_update(agentobj, VB_pos_pd_1D, VB_pos_pd_2D)
        
        % This is a helper function of the agent object's decide method. It
        % evaluates the virtually simulated one-step lookahead observation
        % updated belief state. 
        %
        % Inputs
        %       agentobj:     agent object without virtually simulated 
        %                     observation dependent belief state
        %       VB_pos_pd_1D: 1D virtually simulated position updated 
        %                     belief state - here corresponds to the prior
        %       VB_pos_pd_2D: 2D virtually simulated position updated 
        %                     belief state (empty if one or no remaining 
        %                     target assumed) - here corresponds to the
        %                     prior
        %
        % Outputs
        %       VB_pos_od:    1D virtually simulated observation updated 
        %                     belief state
        %       p_o_vb:       observation probability (evidence)
        %
        % Copyright (C) Lilla Horvath, Dirk Ostwald
        % -----------------------------------------------------------------
            
        % if no target discovery assumed 
        % -----------------------------------------------------------------
        if agentobj.s_t_vb(2) == 0

            % no target found in attempt up to trial
            % -------------------------------------------------------------
            if agentobj.s_t(2) == 0

                % initialize
                p_o_giv_s_vb = squeeze(agentobj.O{1}(agentobj.s_t_vb(1),:,:,:)); % p(o_i = 1|s) (i = 1,2,3,4), NaN for not available observations
                p_s_vb_giv_o = zeros(agentobj.nn,agentobj.nn);              % p(s|[o_1,o_2,o_3,o_4]) initialization

                % cycle over target location pairs
                for l_t1 = 1:agentobj.nn
                    for l_t2 = 1:agentobj.nn

                        % determine emission probabilities based on observation
                        p_oi_giv_s_vb = NaN(1,size(p_o_giv_s_vb,3));

                        % cycle over observations and recover probabilities p(o_i|s)
                        for i = 1:size(p_o_giv_s_vb,3)
                            switch agentobj.o_t_vb(i)
                                case 1
                                    p_oi_giv_s_vb(i) = p_o_giv_s_vb(l_t1,l_t2,i);     % light bar observation
                                case 0                                     
                                    p_oi_giv_s_vb(i) = 1 - p_o_giv_s_vb(l_t1,l_t2,i); % dark bar observation
                            end
                        end

                        % change possible NaN values to 1 for use of prod
                        p_oi_giv_s_vb(isnan(p_oi_giv_s_vb)) = 1;

                        % evaluate Bayes theorem numerator
                        p_s_vb_giv_o(l_t1,l_t2) = prod(p_oi_giv_s_vb)*VB_pos_pd_2D(l_t1,l_t2);
                    end
                end

                % probability of virtually simulated observed bar combination
                p_o_vb = sum(p_s_vb_giv_o(:));

                % normalize
                if p_o_vb ~= 0 
                    p_s_vb_giv_o = p_s_vb_giv_o/sum(p_s_vb_giv_o(:));
                else
                    p_s_vb_giv_o = zeros(agentobj.nn,agentobj.nn);
                end

                % marginalize
                VB_pos_od = sum(p_s_vb_giv_o,2);

            % target found in attempt
            % -------------------------------------------------------------
            else

                % initialize
                p_o_giv_s_vb = squeeze(agentobj.O{2}(agentobj.s_t_vb(1),:,:)); % p(o_i = 1|s) (i = 1,2,3,4), NaN for not available observations
                p_s_vb_giv_o = zeros(agentobj.nn,1);

                % cycle over target locations
                for l_t = 1:agentobj.nn

                    % determine emission probabilities based on observation
                    p_oi_giv_s_vb = NaN(1,size(p_o_giv_s_vb,2));

                    % cycle over observations and recover probabilities p(o_i|s)
                    for i = 1:size(p_o_giv_s_vb,2)
                       switch agentobj.o_t_vb(i)
                            case 1
                                p_oi_giv_s_vb(i) = p_o_giv_s_vb(l_t,i);     % light bar observation
                            case 0
                                p_oi_giv_s_vb(i) = 1 - p_o_giv_s_vb(l_t,i); % dark bar observation
                        end
                    end

                    % change possible NaN values to 1 for use of prod
                    p_oi_giv_s_vb(isnan(p_oi_giv_s_vb)) = 1;

                    % evaluate Bayes theorem numerator
                    p_s_vb_giv_o(l_t) = prod(p_oi_giv_s_vb)*VB_pos_pd_1D(l_t);

                end

                % probability of virtually simulated observed bar combination
                p_o_vb = sum(p_s_vb_giv_o);

                % normalize
                if p_o_vb ~= 0 
                    VB_pos_od = p_s_vb_giv_o/sum(p_s_vb_giv_o);
                else
                    VB_pos_od = zeros(agentobj.nn,1);
                end                    

            end

        % if target discovery assumed 
        % -----------------------------------------------------------------
        elseif agentobj.s_t_vb(2) == 1

            % no target found in attempt up to trial
            % -------------------------------------------------------------
            if agentobj.s_t(2) == 0

                % initialize
                p_o_giv_s_vb = squeeze(agentobj.O{2}(agentobj.s_t_vb(1),:,:)); % p(o_i = 1|s) (i = 1,2,3,4), NaN for not available observations
                p_s_vb_giv_o = zeros(agentobj.nn,1);

                % cycle over target locations
                for l_t = 1:agentobj.nn

                    % determine emission probabilities based on observation
                    p_oi_giv_s_vb = NaN(1,size(p_o_giv_s_vb,2));

                    % cycle over observations and recover probabilities p(o_i|s)
                    for i = 1:size(p_o_giv_s_vb,2)
                       switch agentobj.o_t_vb(i)
                            case 1
                                p_oi_giv_s_vb(i) = p_o_giv_s_vb(l_t,i);     % light bar observation
                            case 0
                                p_oi_giv_s_vb(i) = 1 - p_o_giv_s_vb(l_t,i); % dark bar observation
                        end
                    end

                    % change possible NaN values to 1 for use of prod
                    p_oi_giv_s_vb(isnan(p_oi_giv_s_vb)) = 1;

                    % evaluate Bayes theorem numerator
                    p_s_vb_giv_o(l_t) = prod(p_oi_giv_s_vb)*VB_pos_pd_1D(l_t);

                end

                % probability of virtually simulated observed bar combination
                p_o_vb = sum(p_s_vb_giv_o);

                % normalize
                if p_o_vb ~= 0 
                    VB_pos_od = p_s_vb_giv_o/sum(p_s_vb_giv_o);
                else
                    VB_pos_od = zeros(agentobj.nn,1);
                end
                
            % target found in attempt 
            % -------------------------------------------------------------    
            else
                VB_pos_od = VB_pos_pd_1D;
                p_o_vb    = 0;
            end

        end
        end
    
    end

end

        
        