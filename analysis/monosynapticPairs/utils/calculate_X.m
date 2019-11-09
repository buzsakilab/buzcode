function X = calculate_X(Tlist,hyper_params)
% Calculate covariate bases at the time of presynaptic spike. 
% 
% INPUT
% Tlist = pre (Tlist{1}) and post (Tlist{2}) synaptic spike times (s)
% hyper_params.XTimeSpan = number of seconds with which to divide session
%
% OUTPUT
% X = covariate matrix in basis space (rows = number of presyn spikes, 
%     columns = spline bases)
%
% Author: Abed Ghanbari


% Spline basis over time
mint = min(cellfun(@min,Tlist));
maxt = max(cellfun(@max,Tlist));
% Normalize presynaptic spike times by range of spike times in this pre/post pair
tt = (Tlist{1}-mint)/(maxt-mint);
if ~isfield(hyper_params,'nonstationary_nsplines')
    hyper_params.nonstationary_nsplines = ceil( (maxt-mint)/hyper_params.XTimeSpan );
end

% Splines are centered on presynaptic spike times - every time hyper_params.XTimeSpan worth
% of spikes happens, you put down a spline. The splines are centered where
% they have their peak
Xfr = getCubicBSplineBasis(tt,hyper_params.nonstationary_nsplines,0);
X = Xfr(:,2:end);

