function [synLat,synStrZ,synStrR,sigBounds,synFlip] = FindSynapse(cch,varargin)

% function [synLat,synStrZ,synStrR,sigBounds,synFlip] = FindSynapse(cch,varargin)
%
% FindSynapse - find putatitve monosynaptic connections based on
% cross-correlograms
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


% Adrien Peyrache 2014


% Parameters
bins = 0.5; %in ms
window = 21; %number of bins
alpha = 0.001;
synWin = [0 5]; %window for synaptic connection in ms
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
  end
  
end
%%

b = bins*([0:size(cch,1)-1]-floor(size(cch,1)/2));
%Index of synaptic window
bix = abs(b)<=synWin(2) & abs(b)>=synWin(1);
%Restricted bin vector
bb = b(bix); 


nbExcShifts = round(excSigWin/bins);
if nbExcShifts==0
    error('Inconsistent duration of excitation windows')
end

nbInhShifts = round(inhSigWin/bins);
if nbInhShifts==0
    error('Inconsistent duration of inhibition windows')
end

% Run cch_conv
[dumy, pred, dumy ] = cch_conv(round(cch),window);
predRest = pred(bix,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Putative excitatory synapse


% Estimation of the upper bound of the confidence interval given by maximum
% of the predictor in the synaptic window and evaluated at 1-alpha of the cdf 
% of a poisson distribution
hiBound = poissinv( 1 - alpha, max(pred(bix,:), [], 1 ));
hiBins = bsxfun( @gt, cch( bix, : ), hiBound );

if nbExcShifts>1
    shiftBs = ones(size(hiBins,1)-nbExcShifts+1,size(hiBins,2));
    
    for ii=1:nbExcShifts
        shiftBs = shiftBs & hiBins(ii:end-nbExcShifts+ii,:);
    end
    hiBins = shiftBs;
end

excPairs = find( any(hiBins) );
pkTimes = NaN(length(excPairs),2);
pkStrZ = NaN(length(excPairs),2);
pkStrR = NaN(length(excPairs),2);
flippedE = NaN(length(excPairs),2);

if ~isempty(excPairs)%Difference between actual ccg and upper bound of the confidence
    %interval
    
    diffBins = cch( bix, excPairs);

    % Some ugly code to account for the requirement of consecutive bins
    % to be significant to assess excitatory synapse
    for ii=1:length(excPairs)
        %For this pair, where are the bins which fulfill the requirement?
        signBins = hiBins( :,excPairs(ii) );
        
        % Keep only the first bin in each window
        ds = diff(signBins);
        signBins(2:end) = ds==1;
        signBins = find(signBins);
        
        % Now, we gonna look for peaks at positive or negative latencies (or both)
        % We just need to be sure that there is only one significant bin
        % per sign of latency (this is very unlikely... unless weird
        % parameters were chosen). By default, we'll take the closest from
        % 0. TO BE CHECKED, NOT DEBUGGED
        
        bbs = bb(signBins);
        if sum(bbs<=0)>1
            ix = find(bbs<=0);
            signBins(1:ix(end-1)) = [];
            bbs = bb(signBins);
        end
        if sum(bbs>0)>1
            ix = find(bbs>0);
            signBins(ix(2):end) = [];
        end
        
        % let's go through each of the significant windows ( maximum of 2, one
        % at positive and one at negative latencies)
        for jj=1:length(signBins);
            % In this window, where is the max?
            [pkVal,pkT] = max(diffBins(signBins(jj):signBins(jj)+nbExcShifts-1,ii));
           
            
            % Let's compute the 95th percentile of the expected
            % ccg at this particular bin
            % The strength of connection is then estimated in z-score.

            predPk = predRest(pkT+signBins(jj)-1,excPairs(ii));
            confInt = sqrt(predPk); %Assuming a Poisson process
            strZ = (pkVal-predPk)/confInt;
            strR = pkVal/predPk;
            
            % Associated actual time
            pkT = bb(signBins(jj)+pkT-1);
            
            if pkT>0 %reference cell is presynaptic;
                pkTimes(ii,2) = pkT;
                pkStrZ(ii,2) = strZ;
                pkStrR(ii,2) = strR;
            else %reference cell is postsynaptic
                pkTimes(ii,1) = -pkT;
                pkStrZ(ii,1) = strZ;
                pkStrR(ii,1) = strR;
                flippedE(ii,1) = 1;
            end
        end
    end
else
    pkTimes = [];
    pkStr = [];
end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Putative inhibitory synapse

% Same a.Is above

loBound = poissinv( alpha, min(pred(bix,:), [], 1) );
loBins = bsxfun( @lt, cch( bix, : ), loBound );
if nbInhShifts>1
    shiftBs = ones(size(loBins,1)-nbInhShifts+1,size(loBins,2));
    
    for ii=1:nbInhShifts
        shiftBs = shiftBs & loBins(ii:end-nbInhShifts+ii,:);
    end
    loBins = shiftBs;
end

inhPairs = find( any(loBins) );
trTimes = NaN(length(inhPairs),2);
trStrZ = NaN(length(inhPairs),2);
trStrR = NaN(length(inhPairs),2);
flippedI = NaN(length(excPairs),2);

if ~isempty(inhPairs)
    
    diffBins = cch( bix, inhPairs);
    

    for ii=1:length(inhPairs)
        signBins = loBins( :,inhPairs(ii) );
        
        ds = diff(signBins);
        signBins(2:end) = ds==1;
        signBins = find(signBins);
        
        bbs = bb(signBins);
        if sum(bbs<=0)>1
            ix = find(bbs<=0);
            signBins(1:ix(end-1)) = [];
            bbs = bb(signBins);
        end
        if sum(bbs>0)>1
            ix = find(bbs>0);
            signBins(ix(2):end) = [];
        end
        
        for jj=1:length(signBins);
            % In this window, where is the min?
            [pkVal,pkT] = min(diffBins(signBins(jj):signBins(jj)+nbInhShifts-1,ii));
          
            predPk = predRest(pkT+signBins(jj)-1,inhPairs(ii));
            confInt = sqrt(predPk); %Assuming a Poisson process
            strZ = (pkVal-predPk)/confInt;
            strR = pkVal/predPk;
            pkT = bb(signBins(jj)+pkT-1);
            
            if pkT>0 %reference cell is presynaptic;
                trTimes(ii,2) = pkT;
                trStrZ(ii,2) = strZ;
                trStrR(ii,2) = -strR;
            else %reference cell is postsynaptic
                trTimes(ii,1) = -pkT;
                trStrZ(ii,1) = strZ;
                trStrR(ii,1) = -strR;
                flippedI(ii,1) = 1;
            end
        end
    end

end

%% Output management
synLat = NaN(size(cch,2),2);
synStrZ = NaN(size(cch,2),2);
synStrR = NaN(size(cch,2),2);
synFlip = NaN(size(cch,2),2);

synLat(excPairs,:) = pkTimes;
synStrZ(excPairs,:) = pkStrZ;
synStrR(excPairs,:) = pkStrR;
synFlip(excPairs,:) = flippedE;
% Here, we run through every inhibitory pairs and check that it's not
% already an excitatory connection...

for ii=1:length(inhPairs)
    if ~isnan(trTimes(ii,1))
        if isnan(synLat(inhPairs(ii),1))
            synLat(inhPairs(ii),1) = trTimes(ii,1);
            synStrZ(inhPairs(ii),1) = trStrZ(ii,1);
            synStrR(inhPairs(ii),1) = trStrR(ii,1);
            synFlip(inhPairs,1) = flippedI(ii,1);
        else
            if 0
            figure(1),clf
            plot(b,cch(:,inhPairs(ii)))
            pause
            end
            warning(['Ref. to target cell connection in pair #' num2str(inhPairs(ii)) ' is both excitatory and inhibitory'])
        end
    end
    if ~isnan(trTimes(ii,2))
        if isnan(synLat(inhPairs(ii),2))
            synLat(inhPairs(ii),2) = trTimes(ii,2);
            synStrZ(inhPairs(ii),2) = trStrZ(ii,2);
            synStrR(inhPairs(ii),2) = trStrR(ii,2);
            synFlip(inhPairs,2) = flippedI(ii,2);
        else
            warning(['Target to ref cell connection in pair #' num2str(inhPairs(ii)) ' is both excitatory and inhibitory'])
        end
    end
end

sigBounds = [loBound(:) hiBound(:)];

