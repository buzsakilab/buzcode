function xup = lin_interp(x, upfactor, LOG)
% xup = lin_interp(x, upfactor, LOG)
%
% upfactor is ratio of new/old sampling rates. LOG flag -> interpolate on log scale.
%
% Created by EWS, 2012

if (nargin < 3)
    LOG = 0;
end

if LOG
    x = log10(x);
end

xup = zeros((length(x)-1)*upfactor + 1,1);

orig_inds = 1:upfactor:length(xup);
for i=1:(length(x)-1)
    xup(orig_inds(i)) = x(i);
    xup((orig_inds(i)+1):(orig_inds(i+1)-1)) = x(i) + (1:(upfactor-1))'*(x(i+1)-x(i))/upfactor;
end
xup(end) = x(end);

if (size(x,1) == 1)
    xup = xup';
end

if LOG
    xup = 10.^xup;
end
