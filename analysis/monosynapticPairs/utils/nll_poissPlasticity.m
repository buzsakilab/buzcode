function [f,df,lambda] = nll_poissPlasticity(bta,X,Yy,opts)
% The log of the binned slow rate (opts.rate) is used as a regressor, with the 
% weight fixed to be 1. Only the basis parameters affect the likelihood and its
% gradient.

% Coupling component
cc = X*bta;
% Predicted rate in bins of size dt (can be interpreted as a probability)
lambda = exp(opts.rate + cc) .* opts.dt;
% Negative log likelihood
f = double(nansum(lambda) - nansum(Yy.*log(lambda)));

% Gradient
df = [ nansum(X.*lambda) - sum(X(logical(Yy),:)) ]';

if any(~isfinite(df))
    f = Inf;
end


end