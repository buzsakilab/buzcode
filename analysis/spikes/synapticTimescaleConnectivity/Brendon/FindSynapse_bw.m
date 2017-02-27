function [synLat,synStrZ,synStrR,sigBounds,synStart,synEnd,synFlip] = FindSynapse_bw(cch,varargin)

% function [synLat,synStrZ,synStrR,sigBounds,synStart,synEnd,synFlip] = FindSynapse(cch,varargin)
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
%     Properties    Values
%    -------------------------------------------------------------------------
%     'bins'        size of bins in ms (default = 0.5)
%     'alpha'       p-value (default = 0.001)
%     'synWin'      boundaries of detection window in ms (default = [1 5])
%     'convWin'     size of the convolution window (in number of bins, default = 21)
%     'excSigWin'   duration of actual ccg above threshold for exc. synapses (default = 1ms)
%     'inhSigWin'   duration of actual ccg below threshold for inh. synapses (default = 2ms)
%
% DEPENDENCIES:sort([cstart cend])
%   cch_conv by Eran Stark


% Adrien Peyrache 2014, slight modifications by Brendon Watson 2014


% Default parameters... changed below by function inputs as given
bins = 0.5; %in ms
window = 21; %number of bins
alpha = 0.001;
synWin = [1 4]; %window for synaptic connection in ms
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

b = bins*([0:size(cch,1)-1]-floor(size(cch,1)/2));%timestamps (ms) of bins
%Index of synaptic window
bix = abs(b)<=synWin(2) & abs(b)>=synWin(1);%logical of all bins, 1's where its in the window
%Restricted bin vector
bb = b(bix); %vector of timestamps of only the bins in the synaptic window


nbExcShifts = round(excSigWin/bins);
if nbExcShifts==0
    error('Inconsistent duration of excitation windows')
end

nbInhShifts = round(inhSigWin/bins);
if nbInhShifts==0
    error('Inconsistent duration of inhibition windows')
end

% Run cch_conv
[dummy, pred, dummy ] = cch_conv(round(cch),window);%pred is the prediction based on smoothing (ie based on nearby bins)
predRest = pred(bix,:);%predrest is the preds of the synaptic windows only

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Putative excitatory synapses


% Estimation of the upper bound of the confidence interval given by maximum
% of the predictor in the synaptic window and evaluated at 1-alpha of the cdf 
% of a poisson distribution
hiBound = poissinv( 1 - alpha, max(pred(bix,:), [], 1 ));%high threshold for each pairwise compairson
hiBins = bsxfun( @gt, cch( bix, : ), hiBound );%find bins above the high threshold

if nbExcShifts>1%require multiple consecutive bins to be beyond threshold
    shiftBs = ones(size(hiBins,1)-nbExcShifts+1,size(hiBins,2));
    
    for ii=1:nbExcShifts
        shiftBs = shiftBs & hiBins(ii:end-nbExcShifts+ii,:);
    end
    hiBinsG = shiftBs;%a new version of above-threshold incidents, essentially indicating the first of consecutive above-thresh incidents
end

excPairs = find( any(hiBinsG) );
pkTimes = NaN(length(excPairs),2);
pkStrZ = NaN(length(excPairs),2);
pkStrR = NaN(length(excPairs),2);
pkStart = NaN(length(excPairs),2);
pkEnd = NaN(length(excPairs),2);
flippedE = NaN(length(excPairs),2);

if ~isempty(excPairs)%Difference between actual ccg and upper bound of the confidence
    %interval
    
    evalBins = cch( bix, excPairs);%grabbing bits of CCGs that are from the positive pairs and in the synaptic window

    % Some ugly code to account for the requirement of consecutive bins
    % to be significant to assess excitatory synapse
    for ii=1:length(excPairs)
        
        %For this pair, where are the bins which fulfill the requirement?
        signBins = hiBinsG( :,excPairs(ii) );
        
        % Keep only the first bin in each window... use diff to find first even if multiple windows
        ds = diff(signBins);
        signBins(2:end) = ds==1;
        signBins = find(signBins);
        
        % Now, we gonna look for peaks at positive or negative latencies (or both)
        
        % We just need to be sure that there is only one significant window
        % per sign of latency (this is very unlikely... unless weird
        % parameters were chosen). By default, we'll take the closest from
        % 0. 
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
            [pkVal,pkT] = max(evalBins(signBins(jj):signBins(jj)+nbExcShifts-1,ii));
           
            

            
            
            
            finalbinidx = signBins(jj)+pkT-1;
            predPk = predRest(finalbinidx,excPairs(ii));
            confInt = sqrt(predPk); %Assuming a Poisson process
            strZ = (pkVal-predPk)/confInt;
            strR = pkVal/predPk;
            
            % Associated actual time
            pkT = bb(finalbinidx);
            
            % Grab all above-thresh bins in the same bunch as this one
            dfs = diff(cat(1,0,hiBins( :,excPairs(ii)),0));
                 % firstbin:last bin before the found one that is at the 
                 % leading edge of a contiguous chain (ie so it must be
                 % this chain)
            firstbin = find(dfs(1:end-1)==1 & (1:size(hiBins,1))'<=finalbinidx,1,'last');
%             firstbin = firstbin+1;%correct for the offset from "diff"
                 % lastbin:first bin after the found one that is at the 
                 % trailing edge of a contiguous chain (ie so it must be
                 % this chain)
            lastbin = find(dfs(2:end)==-1 & (1:size(hiBins,1))'>=finalbinidx,1,'first');
            
            
            if pkT>0 %reference cell is presynaptic;
                pkTimes(ii,2) = pkT;
                pkStrZ(ii,2) = strZ;
                pkStrR(ii,2) = strR;
                
                firstbin = max([firstbin find(bb>0,1,'first')]);
                pkStart(ii,2) = bb(firstbin);%make sure start isn't on wrong side of zero
                pkEnd(ii,2) = bb(lastbin);
            else %reference cell is postsynaptic
                pkTimes(ii,1) = -pkT;
                pkStrZ(ii,1) = strZ;
                pkStrR(ii,1) = strR;

                lastbin = min([lastbin find(bb<0,1,'last')]);
                pkStart(ii,1) = bb(lastbin);%note that "FIRST BIN" IS CLOSEST TO ZERO
                pkEnd(ii,1) = bb(firstbin);% "lastbin" is FARTHEST FROM ZERO, so will be first data point from left in flipped cnxns
                flippedE(ii,1) = 1;
            end
        end
    end
else
    pkTimes = [];
    pkStr = [];
    pkStart = [];
    pkEnd = [];
    flippedE = [];
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
    loBinsG = shiftBs;
end

inhPairs = find( any(loBinsG) );
trTimes = NaN(length(inhPairs),2);
trStrZ = NaN(length(inhPairs),2);
trStrR = NaN(length(inhPairs),2);
trStart = NaN(length(inhPairs),2);
trEnd = NaN(length(inhPairs),2);
flippedI = NaN(length(inhPairs),2);

if ~isempty(inhPairs)
    
    evalBins = cch( bix, inhPairs);%grab bins where synapses may have happened
    

    for ii=1:length(inhPairs)
        signBins = loBinsG( :,inhPairs(ii) );
        
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
            [pkVal,pkT] = min(evalBins(signBins(jj):signBins(jj)+nbInhShifts-1,ii));
          
            finalbinidx = signBins(jj)+pkT-1;
            
            predPk = predRest(finalbinidx,inhPairs(ii));
            confInt = sqrt(predPk); %Assuming a Poisson process
            strZ = (pkVal-predPk)/confInt;
            strR = pkVal/predPk;
            
            pkT = bb(finalbinidx);
            
            % Grab all above-thresh bins in the same bunch as this one
            dfs = diff(cat(1,0,loBins( :,inhPairs(ii)),0));
                 % firstbin:last bin before the found one that is at the 
                 % leading edge of a contiguous chain (ie so it must be
                 % this chain)
            firstbin = find(dfs(1:end-1)==1 & (1:size(loBins,1))'<=finalbinidx,1,'last');
%             firstbin = firstbin+1;%correct for the offset from "diff"
                 % lastbin:first bin after the found one that is at the 
                 % trailing edge of a contiguous chain (ie so it must be
                 % this chain)
            lastbin = find(dfs(2:end)==-1 & (1:size(loBins,1))'>=finalbinidx,1,'first');

            
            if pkT>0 %reference cell is presynaptic;
                trTimes(ii,2) = pkT;
                trStrZ(ii,2) = strZ;
                trStrR(ii,2) = -strR;
                
                firstbin = max([firstbin find(bb>0,1,'first')]);
                trStart(ii,2) = bb(firstbin);%make sure start isn't on wrong side of zero
                trEnd(ii,2) = bb(lastbin);
            else %reference cell is postsynaptic
                trTimes(ii,1) = -pkT;
                trStrZ(ii,1) = strZ;
                trStrR(ii,1) = -strR;
                
                lastbin = min([lastbin find(bb<0,1,'last')]);
                trStart(ii,1) = bb(lastbin);
                trEnd(ii,1) = bb(firstbin);
                flippedI(ii,1) = 1;
            end
        end
    end

end

%% Output management
synLat = NaN(size(cch,2),2);
synStrZ = NaN(size(cch,2),2);
synStrR = NaN(size(cch,2),2);
synStart = NaN(size(cch,2),2);
synEnd = NaN(size(cch,2),2);
synFlip = NaN(size(cch,2),2);

if ~isempty(excPairs)
    synLat(excPairs,:) = pkTimes;
    synStrZ(excPairs,:) = pkStrZ;
    synStrR(excPairs,:) = pkStrR;
    synStart(excPairs,:) = pkStart;
    synEnd(excPairs,:) = pkEnd;
    synFlip(excPairs,:) = flippedE;
    % Here, we run through every inhibitory pairs and check that it's not
    % already an excitatory connection...
end

if ~isempty(inhPairs)
    for ii=1:length(inhPairs)
         if ~isnan(trTimes(ii,1))
            if isnan(synLat(inhPairs(ii),1))
                synLat(inhPairs(ii),1) = trTimes(ii,1);
                synStrZ(inhPairs(ii),1) = trStrZ(ii,1);
                synStrR(inhPairs(ii),1) = trStrR(ii,1);
                synStart(inhPairs(ii),1) = trStart(ii,1);
                synEnd(inhPairs(ii),1) = trEnd(ii,1);
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
                synStart(inhPairs(ii),2) = trStart(ii,2);
                synEnd(inhPairs(ii),2) = trEnd(ii,2);
                synFlip(inhPairs,2) = flippedI(ii,2);
            else
                warning(['Target to ref cell connection in pair #' num2str(inhPairs(ii)) ' is both excitatory and inhibitory'])
            end
        end
    end
end

sigBounds = [loBound(:) hiBound(:)];

