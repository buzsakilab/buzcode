function [synLat,synStrZ,synStrR,sigBounds,synStart,synEnd,synFlip] = FindSynapse_SameShank_bw(cch,varargin)

% function [synLat,synStrZ,synStrR,sigBounds,synStart,synEnd,synFlip] = FindSynapse_SameShank_bw(cch,varargin)
%
% FindSynapse_SameShank - find putatitve monosynaptic connections based on
% cross-correlograms for same-shank cell pairs.  Calls FindSynapse.m
% 
% USAGE
% [synLat,synStrZ,synStrR,bounds] = FindSynapse(cch,options)
% 
%     cch: a matrix of cross-correlograms in columns (expressed in number of events)
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
%     sigBounds:   two column matrix of lower (1st column) and upper (2nd
%               column) confidence interval at the given alpha
%     synStart: two column matrix of start times for each found conection
%     synEnd:   two column matrix of start times for each found conection
%     synFlip:  two column binary matrix with 0 in cases where pre was ref
%                  and 1 in cases where the presynaptic cell was target
%   
% <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values(ratio of actual value minus predicted value divided 
%    -------------------------------------------------------------------------
%     'bins'        size of bins in ms (default = 0.5)
%     'alpha'       p-value (default = 0.001)
%     'synWin'      boundaries of detection window in ms (default = [1 5])
%     'convWin'     size of the convolution window (in number of bins, default = 21)
%     'doubleBinSig'  require two consecutive bins to be significant (default = 1)
%     'deadTime'    dead time of spike detection algorithm (default = 0.8ms)
%     'excSigWin'   duration of actual ccg above threshold for exc. synapses (default = 1ms)
%     'inhSigWin'   duration of actual ccg below threshold for inh. synapses (default = 2ms)
%
% DEPENDENCIES:
%   cch_conv by Eran Stark
%   FindSynapse


% Adrien Peyrache 2014


% Parameters

bins = 0.5; %in ms
window = 21; %number of bins
alpha = 0.001;
synWin = [1 4]; %window for synaptic connection in ms
deadTime = 0.8; %in ms
excSigWin = 1; %in ms
inhSigWin = 1.5; %in ms

if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help FindSynapse'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
  if ~isa(varargin{i},'char'),
    error(['Parameter ' num2str(i) ' is not a property (type ''help FindSynapse'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'bins',
      bins = varargin{i+1};
      if ~isa(bins,'numeric') | length(bins) ~= 1 | bins <= 0
           error('Incorrect value for property ''bins'' (type ''help FindSynapse'' for details).');
      end
     case 'alpha',
      alpha = varargin{i+1};
      if ~isa(alpha,'numeric') | length(alpha) ~= 1 | alpha <= 0
          error('Incorrect value for property ''alpha'' (type ''help FindSynapse'' for details).');
      end 
     case 'synwin',
      synWin = varargin{i+1};
      if ~isa(synWin,'numeric') | length(synWin) ~= 2 
           error('Incorrect value for property ''synWin'' (type ''help FindSynapse'' for details).');
      end  
     case 'convwin',
      window = varargin{i+1};
      if ~isa(window,'numeric') | window < 1 
           error('Incorrect value for property ''convWin'' (type ''help FindSynapse'' for details).');
      end
      case 'excsigwin',
      excSigWin = varargin{i+1};
      if ~isa(excSigWin,'numeric')
           error('Incorrect value for property '' excSigWin'' (type ''help FindSynapse'' for details).');
      end
      case 'inhsigwin',
      inhSigWin = varargin{i+1};
      if ~isa(inhSigWin,'numeric')
           error('Incorrect value for property '' inhSigWin'' (type ''help FindSynapse'' for details).');
      end
       case 'deadtime',
      deadtime = varargin{i+1};
      if ~isa(deadtime,'numeric')
           error('Incorrect value for property ''deatime'' (type ''help FindSynapse_SameShank'' for details).');
      end
  end
  
end

b = bins*([0:size(cch,1)-1]-floor(size(cch,1)/2));
%Index of synaptic window
bix = abs(b)<=synWin(2) & abs(b)>=synWin(1);
%Restricted bin vector
bb = b(bix); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correct for dead times in cross-corr

% Large convolution
winLong = 10*window;

if size(cch,1)<2*winLong
    error('CCHs for same shank cell pairs have to be longer than twenty times the convolution window')
end

[dumy, predLong, dumy ] = cch_conv(round(cch),winLong);

% Replace dead bins by prediction over long timescale
binstarts = b-bins/2;
binstops = b+bins/2;
deadBins = find(abs(binstarts)<deadTime & abs(binstops)<deadTime);
cchCor = cch;
cchCor(deadBins,:) = predLong(deadBins,:);

% For partially overlapping bins with detection dead time, renormalize the
% measure
partialBins = find((deadTime>binstarts & deadTime<binstops) | (-deadTime>binstarts & -deadTime<binstops));
correctionRatio = bins/(binstops(partialBins(2))-deadTime);
cchCor(partialBins,:) = round(cchCor(partialBins,:).*repmat(correctionRatio',[length(partialBins) size(cch,2)]));

% Call FindSynapse with the corrected ccgs
[synLat,synStrZ,synStrR,sigBounds,synStart,synEnd,synFlip] = FindSynapse_bw(cchCor,'synWin',synWin,'bins',bins,'alpha',alpha,'excSigWin',excSigWin,'inhSigWin',inhSigWin,'convWin',window);


