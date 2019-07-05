function strength = SpikeTransfer(s1,s2,BinSize,sigbinms)
% this is the core simple measure of spike transfer between a single pair
% of cells... time restriction on S can happen earlier
% s1 is pre, s2 is post
% INPUTS
% s1 = a list of spiketimes for the presynaptic cell
% s2 = a list of spiketimes for the postsynaptic cell
% BinSize = width of bins of interest, in milliseconds
% sigbinms = Signifies the timespan in ms during which connection strength 
%           will be quantified, should be a 2 element vector, element 1 is  
%           the start time element 2 is the end time
% (SampleRate is assumed to be 10000 per TSToolbox)
%
% OUTPUTS
% strength = proportion of presynaptic spikes followed by postsynaptic spikes in the
%            specified time (expressed as a fraction of 1)
%
% Brendon Watson April 2014

% BinSize = 0.5;%in ms
SampleRate = 10000;%%Based on usual TSD object sample rate<< NEED TO GENERALIZE THIS, HAVEN'T FIGURED OUT HOW YET

numSampsBinSize = BinSize*SampleRate/1000;
HalfBins = round(max(abs(sigbinms))*10/numSampsBinSize);

res = [s1' s2'];
clu = [ones(size(s1))' 2*ones(size(s2))'];
%% Get raw CCGs for all pairs (will not use these for same-shank cells
[ccgR, tR] = CCG(res, clu, numSampsBinSize, HalfBins, SampleRate, unique(clu), 'count'); %calc cross correlograms, output as counts... 3D output array

% binstarts = tR-BinSize/2;
% binstops = tR+BinSize/2;
gb = (find(tR==sigbinms(1)):find(tR==sigbinms(2)));

prespikes = length(s1);
%handling various degenerate cases where there were not enough spikes to
%generate a full-dimensioned ccg
if sum(sum(ccgR))==0 
    postspikes = 0;
elseif isempty(ccgR) | isempty(s2);
    postspikes = 0;
elseif prespikes == 0
    postspikes = sum(ccgR(gb));
else
    postspikes = sum(ccgR(gb,1,2));
end

strength = postspikes/prespikes;
if isnan(strength);
    strength = 0;
elseif strength == Inf
    strength = 1;
end

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