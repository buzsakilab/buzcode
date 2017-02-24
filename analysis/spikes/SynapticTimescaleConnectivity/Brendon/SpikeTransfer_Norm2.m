function [strengthbyratio, strengthbyratechg, varargout] = SpikeTransfer_Norm2(times1,times2,binSize,duration,sigbinms,ConvWidth,sameshank,SampleRate)
% Calculates the spiking amplitiude in a pair of cells in bins of interest 
% vs that predicted by the nearby ccg of the same pair.  (Time restriction 
% on S can happen earlier before calling this function).
%
%    s1         first (reference) list of timestamps (in s)
%    s2         optional second list of timestamps (in s)
%    binSize     bin size in s (default = 0.01)
%    duration    duration (in s) of each half of xcorrelogram (default = 2)
%   
% ConvWidth = how widely (#sec) to convolve over when generating estimates of
%               predicted rates (ie 0.006)
% sameshank = binary indicator of whether the given pair on the same shank...
%               determines how the expected rate is calculated (1=same,
%               0=diff)
% % SampleRate = sample rate of timestamps specified in s1 and s2 (default is
% %               10000Hz for TSToolbox)
%
% OUTPUTS
% strengthbyratio = ratio of the bins of interest to the prediciton 
% srengthbyratechg = rate (hz) of extra spikes above/below predicted
% 
% Optional OUTPUT
% ccg - 1D CCG of s1,s2 using specified parameters
%
% Brendon Watson May 2014

% BinSize = 0.5;%in ms

if ~exist('SampleRate','var')
    SampleRate = 10000;%%Based on usual TSD object sample rate<< NEED TO GENERALIZE THIS, HAVEN'T FIGURED OUT HOW YET
end

if isempty(times1) | isempty(times2)
    strengthbyratio = NaN;
    strengthbyratechg = NaN;
else
    times1 = times1(:);
    times2 = times2(:);
    times = [times1;times2];
    groups = 1+[zeros(size(times1));ones(size(times2))];
    HalfBins = round(duration/binSize);
    ConvWidth = ConvWidth/binSize;
%     numSampsBinSize = binSize/1000*SampleRate;%bin size in number of samples
%     HalfBins = duration/1000*SampleRate/numSampsBinSize;
%     % ConvWidth = 0.012; %# of sec
% 
%     res = [s1' s2'];
%     clu = [ones(size(s1))' 2*ones(size(s2))'];

    %% Get raw CCGs for all pairs (will not use these for same-shank cells
    if sameshank% if same shank, do more bins and also smooth wider, but with the same bin width
        HalfBins2 = HalfBins*10;%widen for sameshank
        ConvWidth2 = ConvWidth*10;
        [ccgR, tR] = CCG_In(times,groups,binSize,HalfBins2,SampleRate,[1 2],'hz');
        tR = tR/1000;%convert from ms, for some reason that is default output
        if size(ccgR,1)<2*ConvWidth2
            error('CCGs for same shank cell pairs have to be longer than twenty times the convolution window')
        end
        [dummy, pred, dummy ] = cch_conv(round(ccgR(:,1,2)), ConvWidth2);
        keeperbinidxs = find(tR>=-duration & tR<=duration);
        ccg = ccgR([keeperbinidxs],1,2);

    %     % Replace dead bins by prediction over long timescale
    %     deadBins = find(abs(b)<deadTime);
    % %     cchCor = ccgR;
    %     ccgR(deadBins,:) = pred(deadBins,:);
    % 
    %     % For partially overlapping bins with detection dead time, renormalize the
    %     % measure
    %     partialBins = [deadBins(1)-1 deadBins(end)+1];
    %     correctionRatio = 1./(1-(deadTime-abs(b(partialBins))+bins/2)/bins);
    %     cchCor(partialBins,:) = ccgR(partialBins,:).*repmat(correctionRatio',[1 size(cch,2)]);
    
       
    else
        [ccgR, tR] = CCG_In(times,groups,binSize,HalfBins,SampleRate,[1 2],'hz');
        tR = tR/1000;%convert from ms, for some reason that is default output
        [dummy, pred, dummy ] = cch_conv(round(ccgR(:,1,2)),ConvWidth);
        ccg = ccgR(:,1,2);
    end

    gb = (find(tR==sigbinms(1)/1000):find(tR==sigbinms(2)/1000));

    measured = sum(ccgR(gb,1,2));
    expected = sum(pred(gb));

    strengthbyratio = measured/expected;
    strengthbyratechg = (measured-expected)/length(gb);
    
%% alternative
%     else
%         [ccgR, tR] = CCG(res, clu, numSampsBinSize, HalfBins, SampleRate, unique(clu), 'count'); %calc cross correlograms, output as counts... 3D output array
%         [dummy, pred, dummy ] = cch_conv(round(ccgR(:,1,2)),ConvWidth);
%     end
% 
%     gb = (find(tR==sigbinms(1)):find(tR==sigbinms(2)));
% 
%     measured = sum(ccgR(gb,1,2));
%     expected = sum(pred(gb));
% 
%     % to get hz, divide the summed bin counts by 1) the num of bins, 2)the ref cell spikes (to normalize for number of incidences) 3) the width of the bin in seconds  
%     numbins = length(gb);
%     numprespikes = length(s1);
%     binscalingfactor = (BinSize/1000);
%     strengthbyratechg = (measured-expected)/binscalingfactor/numprespikes/numbins;
end

if nargout==3;
    varargout{1} = ccg;
end




function [out, t, Pairs] = CCG_In(T, G, BinSize, HalfBins, SampleRate, GSubset, Normalization, Epochs)
% Copied from CCG.m from Harris/Sirota...
%
% constructs multiple cross and Auto correlogram
% usage: [ccg, t, pairs] = CCG(T, G, BinSize, HalfBins, SampleRate, GSubset, Normalization, Epochs)
% 
% INPUTS:
% T gives the time of the events, (events need not be sorted by TIME)-in timepoints
% G says which one is in which group
% BinSize gives the size of a bin in input units (i.e. not scaled by SampleRate)
% HalfBins gives the number of bins on each side of 0 - so the total is 1+2*HalfBins
% SampleRate is for x-axis scaling only.  It defaults to 20000
% GSubset says which groups to plot the CCGS of (defaults to all but group 1)
% Normalization indicates the type of y-axis normalization to be used.  
% 'count' indicates that the y axis should show the raw spike count in each bin.
% 'hz' will normalize to give the conditional intensity of cell 2 given that cell 1 fired a spike (default)
% 'hz2' will give the joint intensity, measured in hz^2.
% 'scale' will scale by both firing rates so the asymptotic value is 1.  This gives you
%   the ratio of the coincidence rate to that expected for uncorrelated spike trains
%
% optional input Epochs allows you to compute from only spikes in certain epochs
% and it will bias-correct so you don't see any triangular shape. 
% Warning: if gaps between epochs are shorter than total CCG length, this will mess
% up the edges.
%
% OUTPUTS
% out: 3d array ith the first dim being time lag and the second two 
% specifying the 2 cells in question (in the order of GSubset)
% If there is no output specified, it will plot the CCGs
% t: gives the central timepoint of each bin in milliseconds it seems (*1000)
% pairs: ?indices of pairs at each point of out?
%
% NOTE
% This file calls a C program so your CCG is computed fast.
% to use it, you need to compile mex file from CCGHeart.c
% run : mex -v CCGHeart.c 
% different architectures will produce different extensions of mex file,
% also different versions of matlab link mex file to different libraries
% they are mostly taken into account in the code of CCG.m, but if your
% version or architecture is different from those - modify CCGFun string to
% match the name of the mex file your compiler generated.
% optional output t gives time axis for the bins in ms
% optional output argument pairs gives a nx2 array with the indices of the spikes
% in each train that fall in the CCG.
% written by Ken Harris 
% small editions Anton Sirota
if nargin<5
	SampleRate = 20000;
end
if nargin<6
	GSubset = unique(G);
%	GSubset = setdiff(GSubset, 1);
end
if nargin<7
	Normalization = 'hz';
end;
if nargin<8
    Epochs = [];
end

if length(G)==1
	G = ones(length(T), 1);
	GSubset = 1;
	nGroups = 1;
else
	nGroups = length(GSubset);
end;


% Prepare Res and Clu arrays.
G=G(:);
T=T(:);

if ~isempty(Epochs)
    Included = find(ismember(G,GSubset) & isfinite(T) & WithinRanges(T,Epochs));
    
    % check gaps between epochs are not too short
    GapLen = Epochs(2:end,1) - Epochs(1:(size(Epochs,1)-1),2);
    TooShort = find(GapLen<BinSize*(HalfBins+.5));
    if ~isempty(TooShort)
        fprintf('WARNING: Epochs ');
        fprintf('%d ', TooShort);
        fprintf('are followed by too-short gaps.\n');
    end
 

else
    Included = find(ismember(G,GSubset) & isfinite(T));
    Epochs = [];%[min(T)-1 max(T)+1];
  
end
Res = T(Included);
% if no spikes, return nothing
if length(Res)<=1
    nBins = 1+2*HalfBins;
    out = zeros(nBins, nGroups, nGroups);
    t = 1000*(-HalfBins:HalfBins)*BinSize/SampleRate;
    Pairs = [];
    return
end

% To make the Clu array we need a indexing array, which SUCKS!
G2Clu = full(sparse(GSubset,1,1:nGroups));
Clu = G2Clu(G(Included));


% sort by time
[Res ind] = sort(Res);
Clu = Clu(ind);

% Now call the C program...
ver = version; 
ver = str2num(ver(1:3));
comp = computer;
%check for version and architecture
%LocalName = ['CCGHeart_local.' comp];
%LocalName = ls([matlabroot '/toolbox/local/CCGHeart_local*']);
% LocalName = ls(['/jp03/u12/jagdish/matlab/jp/CCGHeart.mexa64']);
% length(LocalName);
% if LocalName
	CCGFun = 'CCGHeart'; % f$%^ it is the .mex file!~!---jp
% else
%    if ver>=7.2
%        [dd host ] = system('hostname');
%     
%        if strfind(host,'urethane')
%            CCGFun = 'CCGHeart_urethane';
%        else
%            CCGFun = 'CCGHeart7_2';
%        end
%    elseif ver<7.2 & ver>=7 
%             CCGFun = 'CCGHeart7';
%    else
%        CCGFun = 'CCGHeart';
%    end
%    if strcmp(comp,'GLNXA64')  %-------------------- shige commented out
%        CCGFun = [CCGFun '_64'];
%    end                        %----------------------------------------
% end
%fprintf('using %s for version %d\n',CCGFun,ver);
    
    
    
% call the program
nSpikes = length(Res);
if nargout>=3
    % fixing the bug of CCGHeart when no spikes fall withing HalfBins (even
    % for autocorrelogram
    if min(diff(Res))<=BinSize*(HalfBins+1)
        [Counts RawPairs] = feval(CCGFun,Res, uint32(Clu), BinSize, uint32(HalfBins));
        rsRawPairs = reshape(RawPairs, [2 length(RawPairs)/2])';
    else
        warning('pairs cannot be computed - no overlap between spikes withing the range you want');
        rsRawPairs = [];
        Counts = feval(CCGFun,Res, uint32(Clu), BinSize, uint32(HalfBins));
    end
else
    Counts = feval(CCGFun,Res, uint32(Clu), BinSize, uint32(HalfBins));
end
% shape the results
nBins = 1+2*HalfBins;
% if there are no spikes in the top cluster, CCGEngine will produce a output the wrong size
nPresent = max(Clu);
Counts = double(reshape(Counts,[nBins nPresent nPresent]));
if nPresent<nGroups
    % extent array size with zeros
    Counts(nBins, nGroups, nGroups) = 0;
end
    
if nargout>=3
    Pairs = Included(ind(double(rsRawPairs) + 1));
end

% OK so we now have the bin counts.  Now we need to rescale it.

% remove bias due to edge effects - this should be vectorized
if isempty(Epochs)
    Bias = ones(nBins,1);
else
    nTerm = [HalfBins:-1:1 , 0.25 , 1:HalfBins];
	Bias = zeros(nBins,1);
    TotLen = 0;
	for e=1:size(Epochs,1)
        EpochLen = Epochs(e,2)-Epochs(e,1);
        EpochBias = clip(EpochLen - nTerm*BinSize,0,inf)*BinSize;
        Bias = Bias+EpochBias';
        TotLen = TotLen + EpochLen;
	end
    Bias = Bias/TotLen/BinSize;
end

if isempty(Epochs)
      Trange = max(Res) - min(Res); % total time
else
       Trange = sum(diff(Epochs,[],2));
end
t = 1000*(-HalfBins:HalfBins)*BinSize/SampleRate;

% count the number of spikes in each group:
for g=1:nGroups
	nSpikesPerGroup(g) = sum(Clu==g);
end;

% normalize each group
for g1=1:nGroups, for g2=g1:nGroups
	switch Normalization
		case 'hz'
			Factor = SampleRate / (BinSize * nSpikesPerGroup(g1));
			AxisUnit = '(Hz)';
		case 'hz2'
			Factor = SampleRate * SampleRate / (Trange*BinSize);	
			AxisUnit = '(Hz^2)';
		case 'count';
			Factor = 1;
			AxisUnit = '(Spikes)';
		case 'scale'
			Factor = Trange / (BinSize * nSpikesPerGroup(g1) * nSpikesPerGroup(g2));
			AxisUnit = '(Scaled)';
		otherwise
			warning(['Unknown Normalization method ', Normalization]);
	end;
% 	ccg(:,g1,g2) = flipud(Counts(:,g1,g2)) * Factor ./repmat(Bias,[1 nGroups,nGroups]); 
 	ccg(:,g1,g2) = flipud(Counts(:,g1,g2)) * Factor ./Bias; 
	ccg(:,g2,g1) = (Counts(:,g1,g2)) * Factor ./Bias; 
	% now plot, if there is no output argument
	if (nargout==0)
		FigureIndex = g1 + nGroups*(nGroups-g2);
		subplot(nGroups,nGroups,FigureIndex);		
	
		% plot graph
%		bar(1000*(-HalfBins:HalfBins)*BinSize/SampleRate, ccg(:,g1,g2));
		bar(t, ccg(:,g1,g2));

		% label y axis
		if g1==g2
%			ylabel(['ACG ', AxisUnit])	
    		FiringRate = SampleRate * nSpikesPerGroup(g1) / Trange;
%     		Ttitle = sprintf('%d (~%5.2fHz)',GSubset(g1),FiringRate);
% 			title(Ttitle);
    		xlabel('ms');
        else 
        %    set(gca, 'xtick', []);
        end
        if g1==1
            ylabel(sprintf('%d', GSubset(g2)));
%			ylabel(['CCG ', AxisUnit])
% 			Ttitle = sprintf('%d vs %d', GSubset(g1), GSubset(g2));
% 			title(Ttitle);
		end
        if g2==nGroups
    		Ttitle = sprintf('%d (~%5.2fHz)',GSubset(g1),FiringRate);
			title(Ttitle);
        end

		axis tight
	end
end,end;

% only make an output argument if its needed (to prevent command-line spew)
if (nargout>0)
	out = ccg;
end


