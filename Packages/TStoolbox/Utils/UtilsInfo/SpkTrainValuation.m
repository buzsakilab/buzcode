function L = SpkTrainValuation(S,f,t,r)

% L = SpkTrainValuation(S,f,t,r)
% 
% computes log-likelihood of spike train S with intensity function f.
% INPUTS:
%     S: spike timings during test epoch
%     f: a tsd describing the predicted intensity function during test epoch
%     t: time vector of the intensity function (must be same length)
%     r: average firing rate during training epoch
% 
% OUTPUT:
%     L: log-likelihood

% Adrien Peyrache, 2014 (following Harris, 2004)

if length(f) ~= length(t)
    error('f and t must be the same length')
end

dt = median(diff(t));
ft = tsd(t,f/r);
ft = Restrict(ft,S);
d = Data(ft);

intF = sum(f-r).*(dt/10000); %express it in sec!
d = sum(log(d(d>0)));

L = d-intF;