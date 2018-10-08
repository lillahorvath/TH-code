classdef th_taskobj

    % This class definition file specifies the properties and methods of an
    % object that implements the treasure hunt task.
    %
    % Copyright (C) Dirk Ostwald, Lilla Horvath
    % ---------------------------------------------------------------------
    
    % properties of the task object
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        
        nn;                 % number of grid-world nodes
        L1;                 % L1 distances
        A;                  % available action set
        O;                  % observation likelihood function
        sol;                % task solved flag
        att;                % task attempt counter
        t;                  % trial counter
        fnd;                % target found flag
        tkn;                % known target location
        s_t;                % state
        o_t;                % observation state
        a_t;                % agent action state
         
    end
    
    % methods of the task object
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function taskobj = th_taskobj(task_init)

        % This function creates a task object of class th_taskobj based
        % on the initialization input structure.
        %
        %   Inputs
        %       task_init:  structure with task initialization fields
        %
        %   Outputs
        %       taskobj:    task object with initialized values
        %
        % Copyright (C) Dirk Ostwald, Lilla Horvath
        % -----------------------------------------------------------------
        
        % set task properties based on input structure
        taskobj.nn     = task_init.nn;
        taskobj.L1     = task_init.L1;
        taskobj.A      = task_init.A;
        taskobj.O      = task_init.O;
        taskobj.sol    = task_init.sol;
        taskobj.att    = task_init.att;
        taskobj.tkn    = task_init.tkn;
        taskobj.t      = task_init.t;
    
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function taskobj = sample_observation(taskobj)
        
        % This function samples an observation (collection of light and 
        % dark grey bars) of the treasure hunt task based on the task's 
        % observation likelihood function and the current task state.
        %
        %   Inputs
        %       taskobj:    task object with likelihood and current state
        %
        %   Outputs
        %       task:       task object with new observation
        %
        % Copyright (C) Dirk Ostwald, Lilla Horvath
        % -----------------------------------------------------------------
        
        nobs = 4;                                                           % cardinality of the observation set
        
        % select observation likelihood
        if taskobj.s_t(2) == 0
            p_o_giv_s = squeeze(taskobj.O{1}(taskobj.s_t(1),taskobj.s_t(3),taskobj.s_t(4),:));
        else
          
            % first treasure found
            if find(taskobj.fnd) == 1
                p_o_giv_s = squeeze(taskobj.O{2}(taskobj.s_t(1),taskobj.s_t(4),:));
            
            % second treasure found
            elseif find(taskobj.fnd) == 2
                p_o_giv_s = squeeze(taskobj.O{2}(taskobj.s_t(1),taskobj.s_t(3),:));
            end
            
        end
        
        % state dependent observation likelihood
        for i = 1:nobs
            taskobj.o_t(i) = binornd(1,p_o_giv_s(i));                       % sample observations, NaNs for non-sampled observations   
        end
       
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function taskobj = state_update(taskobj, a)
        
        % This function implements the task state update based on the
        % current task state and the current agent action.
        %
        %   Inputs
        %       taskobj:    task object with current state
        %       a:          current agent action \in {-5,+1,+5,-1}
        %
        %   Outputs
        %       taskobj:    task object with udpated state
        %
        % Copyright (C) Dirk Ostwald, Lilla Horvath
        % -----------------------------------------------------------------
        
        % agent action dependent task state update
        taskobj.s_t(1) = taskobj.s_t(1) + a; 

        end
        
    end
end
        
        