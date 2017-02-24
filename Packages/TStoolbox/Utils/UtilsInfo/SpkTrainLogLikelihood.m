function L = SpkTrainLogLikelihood(q,f)

% L = SpkTrainValuation(S,f,r)
% 
% computes log-likelihood of spike train S with intensity function f.
% INPUTS:
%     q: binned spike train
%     f: a tsd describing the predicted intensity function during test epoch
%     dt: time step
% 
% OUTPUT:
%     L: log-likelihood

% Adrien Peyrache, 2014 (following Harris, 2004)

% if length(f) ~= length(q)
%     keyboard
%     error('f and q must be the same length')
% end
% f = f*dt;

L = -f+q.*log(f);
L = sum(L(f>0));