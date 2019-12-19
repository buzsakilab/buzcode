function zerolagrange = CalcZeroLagRange(funcsynapses,pre,post,abovebelow)

% function zerolagrange = CalcZeroLagRange(funcsynapses,pre,post,abovebelow)
%
% CalcZeroLagRange - find start & stop of previously-detected zerolag
% correlation.  Uses bin width/params from previously calc'd zerolag 
%
%  Based on FindSynapse.m
% 
% USAGE
% zerolagrange = CalcZeroLagRange(funcsynapses,pre,post,abovebelow)
% 
%     funcsynapses: a matrix of cross-correlograms in columns (expressed in number of events)
%
% Returns:
%     synLat:   two column vectors of synaptic latencies:
%                 1st column: reference to target cell
%                 2nd column: target to reference cell
%               NaN values indicate non significant interaction
%     synStrZ:  same format as above now indicating synaptic strength in
%               z-scores (ratio of CCG peak minus predicted CCG divided by 
%               root mean square (assuming Poisson dist.) of the predictor)
%               Positive (negative) values indicate excitatory
%               (inhibitory) connections.
%     synStrR:  same as above but synaptic strength is now evaluated as the
%               ratio between the actual CCG peak (or trough) and the predictor
%     bounds:   two column matrix of lower (1st column) and upper (2nd
%               column) confidence interval at the given alpha
%
% <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'bins'        size of bins in ms (default = 0.5)
%     'alpha'       p-value (default = 0.001)
%     'synWin'      boundaries of detection window in ms (default = [1 5])
%     'convWin'     size of the convolution window (in number of bins, default = 21)
%     'excSigWin'   duration of actual ccg above threshold for exc. synapses (default = 1ms)
%     'inhSigWin'   duration of actual ccg below threshold for inh. synapses (default = 2ms)
%
% DEPENDENCIES:
%   cch_conv by Eran Stark


% Brendon Watson 2014


% Parameters
alpha = 0.135;%based on 1.5SD above mean = 68percentile
window = funcsynapses.ZeroLag.ConvolutionWidthInBins;
bins = funcsynapses.ZeroLag.BinMs;

%basic cch/timing info
cch = funcsynapses.ZeroLag.fullCCGMtx(:,pre,post);
b = bins*([0:size(cch,1)-1]-floor(size(cch,1)/2));
zerobin = ceil(size(b,2)/2);

% Run cch_conv
[dumy, pred, dumy ] = cch_conv(round(cch),window);

switch abovebelow
    case 'above'
        Bound = poissinv(1-alpha, max(pred, [], 1));%~1SD... for later measurement of width of any detected bumps
        minnumbins = funcsynapses.ZeroLag.ExcMinimumWidth/bins;
        ab = continuousabove2(cch,Bound,minnumbins,Inf);
    case 'below'
        Bound = poissinv( alpha, min(pred, [], 1));
        minnumbins = funcsynapses.ZeroLag.InhMinimumWidth/bins;
        ab = continuousabove2(cch,Bound,minnumbins,Inf);
end
%find one that spans zero
zerolagrange = ab(ab(:,1)<=zerobin & ab(:,2)>=zerobin,:);
% widthinms = diff(ab)*bins;


% 
% 
% function qualbins = FindConsecBinsIncludingZero(qualBins,abovebelow,nbshifts)
% % spits out bins that are at least the specified duration and which are
% % contiguous with zero.  Not pretty, but does exactly what I want...  Spits
% % out in format the rest of this program is looking for
% 
% for a = 1:size(qualBins,2);
%     data = qualBins(:,a);
%     switch abovebelow
%         case 'above'
%             ab = continuousabove2(data',0.5,1,Inf);
%         case 'below'
%             ab = continuousbelow2(data',0.5,1,Inf);
%     end
%     
%     val = zeros(size(qualBins,1)-1,1);
%     if ~isempty(ab)
%         if bb(ab(1))<=0 && bb(ab(2))>=0
%             if nbshifts>1
%                 shiftBs = ones(size(data,1)-nbshifts+1,size(data,2));
% 
%                 for ii=1:nbshifts
%                     shiftBs = shiftBs & data(ii:end-nbshifts+ii,:);
%                 end
%                 val = shiftBs;
%             end
%         end
%     end
%     qualbins(:,a) = val;
% end
% 
% % hiBins = hiBins2;