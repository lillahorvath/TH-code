function bmselect = th_agent_select(bmselect)

% This function computes the agent model BIC scores based on the model maximum 
% log likelihood values and the number of data points (valid participant 
% % decisions) and subsequently evaluates the hypothesis-averaged posterior 
% model probabilities using SPM's Bayesian model selection routine. 
%
% Inputs
%       bmselect: structure with fields
%        .ll:     agent models maximum log likelihood 
%        .dn:     number of valid participant decisions
%
% Outputs
%       bmselect: input structure with additional fields
%        .bic:    agent model BIC scores
%        .pxp:    hypothesis-averaged posterior model probabilities
% 
% Copyright (C) Lilla Horvath, Dirk Ostwald
% -------------------------------------------------------------------------     

% initialize
% -------------------------------------------------------------------------

% invoke SPM 
addpath(genpath([pwd '\SPM12']))  

% unpack input structure
ll = bmselect.ll;
dn = bmselect.dn;

% perform agent model selection
% -------------------------------------------------------------------------

% ----------------------- compute agent BIC scores ------------------------

% agent with no free parameters: BF-0
bic_BF_0 = 2.*-ll(:,1) + 0.*(log(dn));

% agents with 1 free parameter: BF-R, BB-R-1, BB-R-2, BB-R-25, BB-R-25-2, BB-I
bic_BF_R_to_BB_I = 2.*-ll(:,2:7) +  repmat(1.*(log(dn)),1,6);

% agent with 2 free parameters: BB-H-C
bic_BB_H_C = 2.*-ll(:,8) + 2.*(log(dn));

% agent with 3 free parameters: BB-H-E
bic_BB_H_E = 2.*-ll(:,9) + 3.*(log(dn));

% concatenate
bic = [bic_BF_0, bic_BF_R_to_BB_I, bic_BB_H_C, bic_BB_H_E];

% ---------------- evaluate model exceedance probabilities ----------------

[alpha_all, exp_r_all, xp_all, pxp_all, bor_all] = spm_BMS(-bic);

% set output structure
% -------------------------------------------------------------------------
bmselect.bic = bic;
bmselect.pxp = pxp_all;

end

