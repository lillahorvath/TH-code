function p = th_softmax(v,b)

% This function implements the normalized exponential function for
% mapping action valences onto action probabilities 
%   
% Inputs
%       v: 1 x k array of action valences
%       b: scalar "inverse temperature" parameter
%
% Outputs
%       p: 1 x k array of action probabilities
%
% Copyright (C) Dirk Ostwald, Lilla Horvath
% -------------------------------------------------------------------------

% softmax function evaluation
p = exp(b*v)/sum(exp(b*v));                                                

end