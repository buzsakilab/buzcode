function [synStrZ,synStrR,synStart,synEnd,synWideStart,synWideEnd,sigBounds] = FindZeroLagCorr(cch,varargin)

% function [synLat,synStrZ,synStrR,sigBounds,synFlip] = FindZeroLagCorr(cch,varargin)
%
% FindZeroLagCorr - find zero-lag temporal correlations between pairs of cells.
%  Based on FindSynapse.m
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
synWin = [0 0.5]; %window for synaptic connection in ms
excSigWin = 0.5; %in ms
inhSigWin = 0.5; %in ms

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
% bix = b == synWin;
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

widehiBound = poissinv( 0.865, max(pred(bix,:), [], 1 ));%~1.5SD... for later measurement of width of any detected bumps
widehiBins = bsxfun( @gt, cch( bix, : ), widehiBound );

% 
% for a = 1:size(hiBins,2);%lazy but does exactly what I want...
%     data = hiBins(:,a)';
%     above = continuousabove2(data,0.5,nbExcShifts,Inf);
%     if ~isempty(above)
%         if min(above)<=0 && max(above)>=0
%             1;
%         end
%     end
% end
    
% if nbExcShifts>1
%     shiftBs = ones(size(hiBins,1)-nbExcShifts+1,size(hiBins,2));
%     
%     for ii=1:nbExcShifts
%         shiftBs = shiftBs & hiBins(ii:end-nbExcShifts+ii,:);
%     end
%     hiBins = shiftBs;
% end
% 
hiBinsG = FindConsecBinsIncludingZero(bb,hiBins,'above',nbExcShifts);
widehiBinsG = FindConsecBinsIncludingZero(bb,widehiBins,'above',nbExcShifts);

% %to make sure that zero bin is included in the consecutive bins
% dhi = diff(hiBins,1) %look for zero diffs at zero point
% dhi = dhi>-1;
% zerobinrow = find(bb==0);
% zerobindata = hiBins(zerobinrow,:);
% zdhi = dhi & zerobindata;%find where both zero is hi and diff was zero... ie zero bin was part of a series of 1s
% zx = find(zdhi);
% for a = zx
%     hiBins(:,a) = zeros(size(hiBins,1),1);
% end


% 
% dhi(XX) = 0%set to zeros any column case where it was not
% dhi(xxx) = 0 %set any remaining rows to zero
% 
% %should only have 1's spanning zero now
% sum(dhi) to find good enough columns
% extract times of bins


excPairs = find( any(hiBinsG) );
% excPairs = find(hiBins);
pkTimes = NaN(length(excPairs),2);
pkStrZ = NaN(length(excPairs),2);
pkStrR = NaN(length(excPairs),2);
pkStart = NaN(length(excPairs),2);
pkEnd = NaN(length(excPairs),2);
widepkStart = NaN(length(excPairs),2);
widepkEnd = NaN(length(excPairs),2);
% flippedE = NaN(length(excPairs),2);

if ~isempty(excPairs)%Difference between actual ccg and upper bound of the confidence
    %interval
    
    diffBins = cch( bix, excPairs);

    % Some ugly code to account for the requirement of consecutive bins
    % to be significant to assess excitatory synapse
    for ii=1:length(excPairs)
        %For this pair, where are the bins which fulfill the requirement?
        signBins = hiBinsG( :,excPairs(ii) );
        
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
            
            finalbinidx = signBins(jj)+pkT-1;
            % Associated actual time
%             pkT = bb(finalbinidx);

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
                        
            %same for lowerthresh widebins
            dfs = diff(cat(1,0,widehiBins( :,excPairs(ii)),0));
            wfirstbin = find(dfs(1:end-1)==1 & (1:size(widehiBins,1))'<=finalbinidx,1,'last');
            wlastbin = find(dfs(2:end)==-1 & (1:size(widehiBins,1))'>=finalbinidx,1,'first');

%             
%             if pkT>0 %reference cell is presynaptic;
%                 pkTimes(ii,2) = pkT;
                pkStrZ(ii,2) = strZ;
                pkStrR(ii,2) = strR;
                pkStart(ii,2) = bb(firstbin);
                pkEnd(ii,2) = bb(lastbin);
                widepkStart(ii,2) = bb(wfirstbin);
                widepkEnd(ii,2) = bb(wlastbin);
%             else %reference cell is postsynaptic
%                 pkTimes(ii,1) = -pkT;
%                 pkStrZ(ii,1) = strZ;
%                 pkStrR(ii,1) = strR;
%                 flippedE(ii,1) = 1;
%             end
        end
    end
else
%     pkTimes = [];
    pkStr = [];
end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Putative inhibitory synapse

% Same a.Is above

loBound = poissinv( alpha, min(pred(bix,:), [], 1) );
loBins = bsxfun( @lt, cch( bix, : ), loBound );

wideloBound = poissinv(1-0.865, max(pred(bix,:), [], 1 ));%~1.5SD... for later measurement of width of any detected bumps
wideloBins = bsxfun( @gt, cch( bix, : ), wideloBound );

% if nbInhShifts>1
%     shiftBs = ones(size(loBins,1)-nbInhShifts+1,size(loBins,2));
%     
%     for ii=1:nbInhShifts
%         shiftBs = shiftBs & loBins(ii:end-nbInhShifts+ii,:);
%     end
%     loBins = shiftBs;
% end


loBins = FindConsecBinsIncludingZero(bb,loBins,'below',nbInhShifts);
wideloBins = FindConsecBinsIncludingZero(bb,wideloBins,'below',nbInhShifts);

inhPairs = find( any(loBins) );
trTimes = NaN(length(inhPairs),2);
trStrZ = NaN(length(inhPairs),2);
trStrR = NaN(length(inhPairs),2);
trStart = NaN(length(excPairs),2);
trEnd = NaN(length(excPairs),2);
widetrStart = NaN(length(excPairs),2);
widetrEnd = NaN(length(excPairs),2);
% flippedI = NaN(length(excPairs),2);

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

            finalbinidx = signBins(jj)+pkT-1;
            % Associated actual time
%             pkT = bb(finalbinidx);

            % Grab all above-thresh bins in the same bunch as this one
            dfs = diff(cat(1,0,loBins( :,excPairs(ii)),0));
                 % firstbin:last bin before the found one that is at the 
                 % leading edge of a contiguous chain (ie so it must be
                 % this chain)
            firstbin = find(dfs(1:end-1)==1 & (1:size(loBins,1))'<=finalbinidx,1,'last');
%             firstbin = firstbin+1;%correct for the offset from "diff"
                 % lastbin:first bin after the found one that is at the 
                 % trailing edge of a contiguous chain (ie so it must be
                 % this chain)
            lastbin = find(dfs(2:end)==-1 & (1:size(loBins,1))'>=finalbinidx,1,'first');
                        
            %same for lowerthresh widebins
            dfs = diff(cat(1,0,wideloBins( :,excPairs(ii)),0));
            wfirstbin = find(dfs(1:end-1)==1 & (1:size(wideloBins,1))'<=finalbinidx,1,'last');
            wlastbin = find(dfs(2:end)==-1 & (1:size(wideloBins,1))'>=finalbinidx,1,'first');

%             
%             if pkT>0 %reference cell is presynaptic;
%                 pkTimes(ii,2) = pkT;
                trStrZ(ii,2) = strZ;
                trStrR(ii,2) = strR;
                trStart(ii,2) = bb(firstbin);
                trEnd(ii,2) = bb(lastbin);
                widetrStart(ii,2) = bb(wfirstbin);
                widetrEnd(ii,2) = bb(wlastbin);
%             else %reference cell is postsynaptic
%                 trTimes(ii,1) = -pkT;
%                 trStrZ(ii,1) = strZ;
%                 trStrR(ii,1) = -strR;
%                 flippedI(ii,1) = 1;
%             end
        end
    end

end

%% Output management
% synLat = NaN(size(cch,2),2);
synStrZ = NaN(size(cch,2),2);
synStrR = NaN(size(cch,2),2);
synStart = NaN(size(cch,2),2);
synEnd = NaN(size(cch,2),2);
synWideStart = NaN(size(cch,2),2);
synWideEnd = NaN(size(cch,2),2);
% synFlip = NaN(size(cch,2),2);

% synLat(excPairs,:) = pkTimes;
synStrZ(excPairs,:) = pkStrZ;
synStrR(excPairs,:) = pkStrR;
synStart(excPairs,:) = pkStart;
synEnd(excPairs,:) = pkEnd;
synWideStart(excPairs,:) = widepkStart;
synWideEnd(excPairs,:) = widepkEnd;
% synFlip(excPairs,:) = flippedE;
% Here, we run through every inhibitory pairs and check that it's not
% already an excitatory connection...

for ii=1:length(inhPairs)
    if ~isnan(trTimes(ii,1))
        if isnan(synLat(inhPairs(ii),1))
%             synLat(inhPairs(ii),1) = trTimes(ii,1);
            synStrZ(inhPairs(ii),1) = trStrZ(ii,1);
            synStrR(inhPairs(ii),1) = trStrR(ii,1);
            synStart(inhPairs,:) = trStart(ii,1);
            synEnd(inhPairs,:) = trEnd(ii,1);
            synWideStart(inhPairs,:) = widetrStart(ii,1);
            synWideEnd(inhPairs,:) = widetrEnd(ii,1);
%             synFlip(inhPairs,1) = flippedI(ii,1);
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
%             synLat(inhPairs(ii),2) = trTimes(ii,2);
            synStrZ(inhPairs(ii),2) = trStrZ(ii,2);
            synStrR(inhPairs(ii),2) = trStrR(ii,2);
            synStart(inhPairs,:) = trStart(ii,2);
            synEnd(inhPairs,:) = trEnd(ii,2);
            synWideStart(inhPairs,:) = widetrStart(ii,2);
            synWideEnd(inhPairs,:) = widetrEnd(ii,2);
%             synFlip(inhPairs,2) = flippedI(ii,2);
        else
            warning(['Target to ref cell connection in pair #' num2str(inhPairs(ii)) ' is both excitatory and inhibitory'])
        end
    end
end

sigBounds = [loBound(:) hiBound(:)];





function qualbins = FindConsecBinsIncludingZero(bb,qualBins,abovebelow,nbshifts)
% spits out bins that are at least the specified duration and which are
% contiguous with zero.  Not pretty, but does exactly what I want...  Spits
% out in format the rest of this program is looking for
% Custom written to be in this code, not supposed to be modular 
%
% INPUTS:
% qualBins = is the input data and is a binary matrix of above-thresh bins
% bb       = is a vector, same length as qualBins, timestamp for each bin
% abovebelow = characterstring which should either be "above" or "below" 
%               specifying whether this is looking for high or low values
% nbshifts = number of bins needed to be significant, from adrien's code
%
% OUTPUT:
% qualbins = logical matrix with 1s at bins above thresh and contiguous
%               other bins of nbshifts and contig with zero time in bb


for a = 1:size(qualBins,2);
    data = qualBins(:,a);
    switch abovebelow
        case 'above'
            ab = continuousabove2(data',0.5,nbshifts,Inf);
        case 'below'
            ab = continuousbelow2(data',0.5,nbshifts,Inf);
    end
    
    val = zeros(size(qualBins,1)-1,1);
    if ~isempty(ab)
        if bb(ab(1))<=0 && bb(ab(2))>=0
            if nbshifts>1
                shiftBs = ones(size(data,1)-nbshifts+1,size(data,2));

                for ii=1:nbshifts
                    shiftBs = shiftBs & data(ii:end-nbshifts+ii,:);
                end
                val = shiftBs;
            end
        end
    end
    qualbins(:,a) = val;
end

% hiBins = hiBins2;