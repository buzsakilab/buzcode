function [strengthbyratio, strengthbyratechg, varargout] = SpikeTransfer_Norm2(s1,s2,BinMs,HalfCCGWidthinms,sigbinms,ConvWidth,sameshank,SampleRate)
% Calculates the spiking amplitiude in a pair of cells in bins of interest 
% vs that predicted by the nearby ccg of the same pair.  (Time restriction 
% on S can happen earlier before calling this function).
% 
% INPUTS
% s1 = a list of spiketimes for the presynaptic cell (in sample num)... can
%             be in seconds if SampleRate = 1;
% s2 = a list of spiketimes for the postsynaptic cell (in sample num)... can
%             be in seconds if SampleRate = 1;
% BinMs = width of bins of interest, in milliseconds
% HalfCCGWidthinms = Number of bins to calculate on one side of the ccg...
%               ccg calculated as part of this script to allow quantification
% sigbinms = Signifies the timespan in ms during which connection strength 
%           will be quantified, should be a 2 element vector, element 1 is  
%           the start time element 2 is the end time
% ConvWidth = how widely (#bins) to convolve over when generating estimates of
%               predicted rates (ie 12)
% sameshank = binary indicator of whether the given pair on the same shank...
%               determines how the expected rate is calculated (1=same,
%               0=diff)
% SampleRate = sample rate of timestamps specified in s1 and s2 (default is
%               10000Hz for TSToolbox)... if SampleRate = 1 then s1 and s2
%               can be in seconds
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

if isempty(s1) | isempty(s2)
    strengthbyratio = NaN;
    strengthbyratechg = NaN;
else
    numSampsBinSize = BinMs/1000*SampleRate;%bin size in number of samples
    HalfBins = HalfCCGWidthinms/1000*SampleRate/numSampsBinSize;
    % ConvWidth = 12; %# of bins

    res = [s1' s2'];
    clu = [ones(size(s1))' 2*ones(size(s2))'];

    %% Get raw CCGs for all pairs (will not use these for same-shank cells
    if sameshank% if same shank, do more bins and also smooth wider, but with the same bin width
        HalfBins2 = HalfBins*10;%widen for sameshank
        ConvWidth2 = ConvWidth*10;
        [ccgR, tR] = CCG(res, clu, numSampsBinSize, HalfBins2, SampleRate, unique(clu), 'hz'); %calc cross correlograms, output as counts... 3D output array
        if size(ccgR,1)<2*ConvWidth2
            error('CCGs for same shank cell pairs have to be longer than twenty times the convolution window')
        end
        [dummy, pred, dummy ] = cch_conv(round(ccgR(:,1,2)), ConvWidth2);
        keeperbinidxs = find(tR>=-HalfCCGWidthinms & tR<=HalfCCGWidthinms);
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
        [ccgR, tR] = CCG(res, clu, numSampsBinSize, HalfBins, SampleRate, unique(clu), 'hz'); %calc cross correlograms, output as counts... 3D output array
        [dummy, pred, dummy ] = cch_conv(round(ccgR(:,1,2)),ConvWidth);
        ccg = ccgR(:,1,2);
    end

    gb = (find(tR==sigbinms(1)):find(tR==sigbinms(2)));

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




% % binstarts = tR-BinSize/2;
% % binstops = tR+BinSize/2;
% 
% prespikes = length(s1);
% %handling various degenerate cases where there were not enough spikes to
% %generate a full-dimensioned ccg
% if sum(sum(ccgR))==0 
%     postspikes = 0;
% elseif isempty(ccgR) | isempty(s2);
%     postspikes = 0;
% elseif prespikes == 0
%     postspikes = sum(ccgR(gb));
% else
%     postspikes = sum(ccgR(gb,1,2));
% end
% 
% strength = postspikes/prespikes;
% if isnan(strength);
%     strength = 0;
% elseif strength == Inf
%     strength = 1;
% end

% 
% 
% % Loop over found synapses... maybe over classes... ie nl, wide, zerolag, 
% sigbins = signif bins in each cnxn
% pre = 
% post = 
% pretrain = Data(TimePoints(S{pre});
% posttrain = Data(TimePoints(S{post});
% numprespikes = length(pre);
% 
% a) find spikes of pre
%    Simply count spikes from train 1
% 
%     b) find post spikes at at significant bins...
%         just take funcsyn.fullCCGmtx
% 
%     %% >> when want to change so SpikeTransferVsTime.m use below
% %     startstops = [min([pretrain(1) postrain(1)] max([pretrain(1) postrain(1)]];%for now just do whole thing
% %     b: [ccg,x,y] = SpecifiedTimeCCG_bw(pretrain,posttrain,startstops,'mode','count','duration',0.03)
% %         ... 'mode','count' for raw number
% %         .... 'duration',0.03 for 30ms
% %         >> take output, keep desired only, bins, add counts
% 
% c: divide b by a
% end