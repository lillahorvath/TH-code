function a_fit = th_agent_fit(a_fit)

% This function implements a maximum log likelihood approach to the
% evaluation of the agent-data log likelihood.
%
% Inputs
%       a_fit:         structure with fields
%        .pomdp:       pomdp task formulation
%        .y:           observable data (observed state-observation-action evolution)
%        .nruns:       number of runs
%        .ntasks:      number of tasks per run
%        .a_class:     agent class
%        .beta_int:    beta interval endpoints 
%        .beta_init:   beta initial point 
%        .lambda_int:  lambda interval endpoints 
%        .lambda_init: lambda initial point 
%        .os_int:      decay offset interval endpoints 
%        .os_init      decay offset initial point 
%        .sl_int:      decay slope interval endpoints 
%        .sl_init      decay slope initial point 
%        .sdir:        participant behavioral results (ll figures) directory
%        .plot_ll:     log likelihood function plot flag
%
% Outputs
%       a_fit:         input structure with additional fields
%        .ll:          maximum log likelihood
%        .n_dp:        number of valid data points
%        .beta_ml:     ML beta estimate
%        .lambda_ml:   ML lambda estimate
%        .os_ml:       ML decay offset estimate
%        .sl_ml:       ML decay slope estimate
%
% Copyright (C) Dirk Ostwald, Lilla Horvath
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%                              initialization 
% -------------------------------------------------------------------------

% transfer input structure
a_ll = a_fit;

% -------------------------------------------------------------------------
%                              evaluate agents 
% -------------------------------------------------------------------------

% ************ direct log likelihood evaluation for agent BF-0 ************

if isequal(a_ll.a_class, @th_agentobj_BF_0)

    % evaluate log likelihood without parameter
    % --------------------------------------------------------------------- 
    
    % set input values
    par      = [];                                                          % no parameters
    mll_flag = 1;                                                           % mll_flag on: return number of valid data points
    
    % evaluate and set log likelihood and number of valid data points
    [ll, n_dp] = eval_agent_ll(a_ll, par, mll_flag);
    a_fit.ll   = ll;
    a_fit.n_dp = n_dp;

% ******** ML-based log likelihood evaluation for all other agents ********

else
    
    % obtain ML estimate: apply constrained nonlinear optimization routine
    % ---------------------------------------------------------------------
    
    % initialize optimization for beta 
    beta_min  = a_ll.beta_int(1);                                           % lower interval endpoint        
    beta_max  = a_ll.beta_int(2);                                           % upper interval endpoint
    beta_init = a_ll.beta_init;                                             % initial point
    
    % set mll_flag to 0 - don't return number of the valid data points
    mll_flag = 0;
    
    % BB-H-C
    if isequal(a_ll.a_class, @th_agentobj_BB_H_C)  
        
        % initialize optimization for lambda 
        lambda_min  = a_ll.lambda_int(1);                                   % lower interval endpoint      
        lambda_max  = a_ll.lambda_int(2);                                   % upper interval endpoint
        lambda_init = a_ll.lambda_init;                                     % initial point
        
        % run optimization routine (beta, lambda)
        par_ml = fmincon(@(par) -eval_agent_ll(a_ll, par, mll_flag), [beta_init, lambda_init], [], [],[], [], [beta_min, lambda_min], [beta_max, lambda_max]);
    
    % BB-H-E
    elseif isequal(a_ll.a_class, @th_agentobj_BB_H_E)  
        
        % initialize optimization for offset and slope 
        os_min  = a_ll.os_int(1);                                           % lower interval endpoint  
        os_max  = a_ll.os_int(2);                                           % upper interval endpoint  
        os_init = a_ll.os_init;                                             % initial point  
        sl_min  = a_ll.sl_int(1);                                           % lower interval endpoint  
        sl_max  = a_ll.sl_int(2);                                           % upper interval endpoint  
        sl_init = a_ll.sl_init;                                             % initial point  
        
        % run optimization routine (beta, offset, slope)
        par_ml = fmincon(@(par) -eval_agent_ll(a_ll, par, mll_flag), [beta_init, os_init, sl_init], [], [],[], [], [beta_min, os_min, sl_min], [beta_max, os_max, sl_max]);    
    
    % all other agents
    else
        
        % run optimization routine (beta)
        par_ml = fmincon(@(par) -eval_agent_ll(a_ll, par, mll_flag), beta_init, [], [],[], [], beta_min, beta_max);
    end
    
    % evaluate maximum log likelihood
    % ---------------------------------------------------------------------
    
    % set input values
    a_fit.beta_ml = par_ml(1);                                              % ML beta estimate
    mll_flag      = 1;                                                      % mll_flag on: return number of valid data points                                                           
    
    % BB-H-C
    if isequal(a_ll.a_class, @th_agentobj_BB_H_C)
        
        % set ML lambda estimate
        a_fit.lambda_ml = par_ml(2);
        
        % evaluate and set log likelihood and number of valid data points
        [ll, n_dp] = eval_agent_ll(a_ll, [a_fit.beta_ml, a_fit.lambda_ml], mll_flag);
        a_fit.ll   = ll;
        a_fit.n_dp = n_dp; 
    
    % BB-H-E    
    elseif isequal(a_ll.a_class, @th_agentobj_BB_H_E)  
        
        % set ML offset and slope estimates
        a_fit.os_ml = par_ml(2);
        a_fit.sl_ml = par_ml(3);
        
        % evaluate and set log likelihood and number of valid data points
        [ll, n_dp] = eval_agent_ll(a_ll, [a_fit.beta_ml, a_fit.os_ml, a_fit.sl_ml], mll_flag);
        a_fit.ll   = ll;
        a_fit.n_dp = n_dp; 
    
    % all other agents    
    else 
        
        % evaluate and set log likelihood and number of valid data points
        [ll, n_dp] = eval_agent_ll(a_ll, a_fit.beta_ml, mll_flag);
        a_fit.ll   = ll;
        a_fit.n_dp = n_dp;
    end
    
    % visualization
    % ---------------------------------------------------------------------
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % TO DO: IF WE WANT THESE FIGS IN THE PAPER (SUPPL), WE NEED TO DO:
    % 1. create Dirk's print version 
    % 2. create Lilla's export_fig version
    % 3. choose 1 or 2
    % 4. visualize BB-H-E (3D)
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if a_fit.plot_ll
        
        % create subject directory containing the figures
        if exist(a_fit.sdir, 'dir')                                         % delete if pre-existent and create a new one
            rmdir(a_fit.sdir,'s')
            mkdir(a_fit.sdir)
        else
            mkdir(a_fit.sdir)                                               % create if non-existent
        end
        
        % initialize log likelihood function support (beta)
        beta_res = 20;
        beta     = linspace(beta_min,beta_max,beta_res);
        ll_beta  = NaN(beta_res,1);
        
        % turn off mll_flag
        mll_flag = 0;
        
        % BB-H-C
        if isequal(a_ll.a_class, @th_agentobj_BB_H_C) 
            
            % initialize log likelihood function support (beta and lambda)
            lambda_res     = 20;
            lambda         = linspace(lambda_min,lambda_max,lambda_res); 
            ll_beta_lambda = NaN(beta_res, lambda_res);
            
            % cycle over beta parameter values
            for b = 1:beta_res
                
                % cycle over lambda parameter values
                for l = 1:lambda_res
                    
                    % evaluate the log likelihood function
                    ll_beta_lambda(b,l) = eval_agent_ll(a_ll, [beta(b),lambda(l)], mll_flag);
                end
            end
            
            % visualize log likelihood matrix
            h = figure;
            set(h, 'Color', [1 1 1])
            a_str = 'BB-H-C';
            imagesc(ll_beta_lambda')
            title(a_str, 'FontName', 'Arial', 'FontSize', 6)
            xlabel('$\beta$', 'Interpreter', 'Latex', 'FontSize', 6)
            ylabel('$\lambda$', 'Interpreter', 'Latex', 'FontSize', 6)
            set(gca, 'xtick', 0:2:beta_res, 'ytick', 0:2:lambda_res, 'xticklabel', round(beta(1:2:end),2), 'yticklabel', round(lambda(1:2:end),2))
            c = colorbar;
            ylabel(c, '$ln(p_{\beta,\lambda}(y))$', 'Interpreter', 'Latex', 'FontSize', 6)
            grid on
            saveas(h, fullfile(a_fit.sdir, [a_str '.fig']))
            close(h)           
        
        % BB-H-E
        elseif isequal(a_ll.a_class, @th_agentobj_BB_H_E) 
            
            return % CHANGE THIS: visualize BB-H-E (3D)!!
        
        % all other agents    
        else

            % over beta parameter values
            for b = 1:beta_res
                
                % evaluate the log likelihood function
                ll_beta(b) = eval_agent_ll(a_ll, beta(b), mll_flag);
            end

            % visualize log likelihood matrix
            h = figure;
            set(h, 'Color', [1 1 1])
            if isequal(a_ll.a_class, @th_agentobj_BF_0)
                a_str = 'BF-0';
            elseif isequal(a_ll.a_class, @th_agentobj_BF_R)
                a_str = 'BF-R';
            elseif isequal(a_ll.a_class, @th_agentobj_BB_R_1)
                a_str = 'BB-R-1';
            elseif isequal(a_ll.a_class, @th_agentobj_BB_R_2)
                a_str = 'BB-R-2';
            elseif isequal(a_ll.a_class, @th_agentobj_BB_R_25)
                a_str = 'BB-R-25';
            elseif isequal(a_ll.a_class, @th_agentobj_BB_R_25_2)
                a_str = 'BB-R-25-2';
            elseif isequal(a_ll.a_class, @th_agentobj_BB_I)
                a_str = 'BB-I';
            end
            plot(beta,ll_beta, 'LineWidth', 1)
            plot(a_fit.beta_ml, a_fit.ll, 'bo', 'MarkerFaceColor', 'b')
            title(a_str, 'FontName', 'Arial', 'FontSize', 6)
            xlabel('$\beta$', 'Interpreter', 'Latex', 'FontSize', 6)
            ylabel('$ln(p_{\beta}(y))$', 'Interpreter', 'Latex', 'FontSize', 6)
            grid on
            saveas(h,fullfile(a_fit.sdir, [a_str '.fig']))
            close(h)
            
        end
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % END OF TO DO
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end
end

function [ll, sum_dp] = eval_agent_ll(a_ll, par, mll_flag)

% This function evaluates the log likelihood of the choice data for a given
% agent based on the observed data history (state-observation-action 
% evolution) and a given set of parameter values
%
% Inputs
%       a_ll:      structure with fields
%        .pomdp:       pomdp task formulation
%        .y:           observable data (observed state-observation-action evolution)
%        .nruns:       number of runs
%        .ntasks:      number of tasks per run
%        .a_class:     agent class
%       par:       set of parameter values (BF-0: empty, BB-H-C: beta, lambda, BB-H-C: beta, offset, slope, all other agents: beta)
%       mll_flag:  evaluate maximum log likelihood flag
%
%  Outputs
%       ll:        choice data log likelihood
%       sum_dp:    number of valid data points
%
% Copyright (C) Dirk Ostwald, Lilla Horvath
% -------------------------------------------------------------------------

% initialization
% -------------------------------------------------------------------------

% unpack input structure
nruns  = a_ll.nruns;
ntasks = a_ll.ntasks;

% initialize log likelihood and data points counter arrays
ll_agent = NaN(nruns, ntasks);
num_dp   = NaN(nruns, ntasks);

% agent evaluation with a given set of parameter values
% -------------------------------------------------------------------------

% cycle over runs
for r = 1:nruns

    % cycle over tasks 
    for t = 1:ntasks

        % specify evaluation input structure
        a_eva         = [];                                                 % initialize 
        a_eva.pomdp   = a_ll.pomdp;                                         % pomdp formulation
        a_eva.y       = a_ll.y{r}{t};                                       % observable data (observed state-observation-action evolution)
        a_eva.a_class = a_ll.a_class;                                       % agent class
        a_eva.par     = par;                                                % agent parameter set
        
        % evaluate agent
        a_eva = th_agent_eva(a_eva); 
        
        % set agent's choice data log likelihood
        ll_agent(r,t) = a_eva.ll;
        
        % set number of valid data points (valid actions)
        num_dp(r,t) = a_eva.nva;
    end
end
    
% evaluate and set output arguments
% -------------------------------------------------------------------------

% sum log likelihoods
ll = sum(ll_agent(:));

% return number of valid data points for ll-based BIC evaluation 
if mll_flag
    sum_dp = sum(num_dp);                                                  
else
    sum_dp = [];
end

end



