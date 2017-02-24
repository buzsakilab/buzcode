function [funcsynapses] = Make_FindZeroLagCorr_bw(S,shank,funcsynapses)

%% Parameters
alpha = 0.01;%signifiance cutoff/pvalue
BinSize = 0.25;%in ms
synapticwindow = [0 .25];%Milliseconds in which synapses much occur (inclusive)...
                        % = synWin input to FindSynapses.m
convolutionwidth = round(6 * 1/BinSize); % size of the convolution window (in number of bins
                        % = convWin input to FindSynapses.m
excSigWin = .5; % duration in ms of actual ccg above threshold for exc. synapses
inhSigWin = .5; % duration in ms of actual ccg below threshold for inh. synapses

[T,G] = oneSeries(S);
T = TimePoints(T);

SampleRate = 10000;%%Based on usual TSD object sample rate<< NEED TO GENERALIZE THIS, HAVEN'T FIGURED OUT HOW YET
numBinsBinSize = BinSize*SampleRate/1000;
HalfBins = round(300/numBinsBinSize);


%% Get raw CCGs for all pairs (will not use these for same-shank cells
log = 0;
if exist('funcsynapses','var')
    if  isfield (funcsynapses,'ZeroLag') 
        if ~isempty(funcsynapses.ZeroLag.fullCCGMtx);
            log = 1;
        end
    end
end
if log == 1
    ccgR = funcsynapses.ZeroLag.fullCCGMtx;
    tR = funcsynapses.ZeroLag.CCGbins;
else
    [ccgR, tR, Pairs] = CCG(T, G, numBinsBinSize, HalfBins, SampleRate, unique(G), 'count'); %calc cross correlograms, output as counts... 3D output array
end

%% Iterate through all pairs
cch = [];
cchSame = [];

cellPairId = [];
cellPairIdSame = [];

for x=1:length(S)
%     rgx = Range(S{x});
%     rx = length(Range(S{x}));
    for y=x+1:length(S)
%        rgy = Range(S{y});

       if shank(x)~=shank(y)
%            [h,b] = CrossCorr(rgx,rgy,bin,nbBins);
%            cch = [cch h*rx*bin/1000];[funcsynapses] = Make_FindSynapse_bw(S,shank)
           cch(:,end+1) = ccgR(:,x,y);
           cellPairId = [cellPairId;[x y]];
%        else
% %            [h,b] = CrossCorr(rgx,rgy,bin,nbBinsLg);
% %            cchSame = [cchSame h*rx*bin/1000];
%            T = cat(1,TimePoints(S{x}),TimePoints(S{y}));
%            G = cat(1,ones(length(S{x}),1),2*ones(length(S{y}),1));
%             
%            [ccgR_temp, tR_temp, Pairs_temp] = CCG(T, G, numBinsBinSize, HalfBins*10, SampleRate, unique(G), 'count'); %calc cross correlograms, output as counts... 3D output array
%            cchSame(:,end+1) = ccgR_temp(:,1,2);
%            cellPairIdSame = [cellPairIdSame;[x y]];
       end       
    end
end 

[synStrZ,synStrR,synStart,synEnd,synWideStart,synWideEnd,bounds] = ...
    FindZeroLagCorr(cch,'bins',BinSize,...
    'synWin',synapticwindow,'alpha',alpha,'convWin',convolutionwidth,...
    'excSigWin',excSigWin,'inhSigWin',inhSigWin);
synWeightsZ = pair2mat(synStrZ,cellPairId,length(S));
synWeightsR = pair2mat(synStrR,cellPairId,length(S));
synSig = pair2mat(bounds,cellPairId,length(S));
synStarts = pair2mat(synStart,cellPairId,length(S));
synEnds = pair2mat(synEnd,cellPairId,length(S));
synWideStarts = pair2mat(synWideStart,cellPairId,length(S));
synWideEnds = pair2mat(synWideEnd,cellPairId,length(S));

% 
% [synLat,synStrZ,synStrR,bounds,synFlip] = FindSynapse_SameShank(cchSame,'bins',BinSize,'alpha',0.01,'synwin',synapticwindow,'convWin',convolutionwidth);
% synTimes = pair2mat(synLat,cellPairIdSame,synTimes);
% synWeightsZ = pair2mat(synStrZ,cellPairIdSame,synWeightsZ);
% synWeightsR = pair2mat(synStrR,cellPairIdSame,synWeightsR);
% synSig = pair2mat(bounds,cellPairIdSame,synSig);
% synRev = pair2mat(synFlip,cellPairIdSame,synRev);
% synRev(isnan(synRev)) = 0;

%% Flipping outputs of above so pre-synaptic cells are on dim 1, Post-syn on dim 2
% synTimes = synTimes';
synWeightsZ = synWeightsZ';
synWeightsR = synWeightsR';
synSig = synSig';
synStarts = synStarts';
synEnds = synEnds';
synWideStarts = synWideStarts';
synWideEnds = synWideEnds';
% synRev = synRev';

%% get total connection strengths based on start and end times of each interaction
% [pre,post] = find(~isnan(synWeightsZ));
% strengths = spikeTransfer_InPairSeries(S,[pre post],BinSize,synStarts,synEnds);
% synTransfers = nan(size(synStarts));
% for a = 1:length(pre)
%     synTransfers(pre(a),post(a)) = strengths(a);
% end

%% Prep for output
funcsynapses.ZeroLag.BinMs = BinSize;
funcsynapses.ZeroLag.AlphaValue = alpha;
% funcsynapses.ZeroLag.SynapticWindowInMs = synapticwindow;
funcsynapses.ZeroLag.ConvolutionWidthInBins = convolutionwidth;
funcsynapses.ZeroLag.ExcMinimumWidth = excSigWin;
funcsynapses.ZeroLag.InhMinimumWidth = inhSigWin;
funcsynapses.ZeroLag.CellShanks = shank;
% funcsynapses.ZeroLag.CellShankIDs = cellIx;
% funcsynapses.FlippedCnxns = synRev;

% funcsynapses.CnxnTimesVsRefSpk = synTimes;
% funcsynapses.CnxnBinsVsRefSpk = funcsynapses.CnxnTimesVsRefSpk./funcsynapses.BinMs;
% funcsynapses.CnxnTimesVsCCGStart = synTimes+funcsynapses.BinMs*floor(length(tR)/2);
% funcsynapses.CnxnBinsVsCCGStart = ceil(length(tR)/2)+funcsynapses.CnxnTimesVsRefSpk./funcsynapses.BinMs;
funcsynapses.ZeroLag.CnxnStartTimesVsRefSpk = synStarts;
funcsynapses.ZeroLag.CnxnStartBinsVsRefSpk = funcsynapses.ZeroLag.CnxnStartTimesVsRefSpk./funcsynapses.ZeroLag.BinMs;
funcsynapses.ZeroLag.CnxnStartTimesVsCCGStart = synStarts+funcsynapses.ZeroLag.BinMs*floor(length(tR)/2);
funcsynapses.ZeroLag.CnxnStartBinsVsCCGStart = ceil(length(tR)/2)+funcsynapses.ZeroLag.CnxnStartTimesVsRefSpk./funcsynapses.ZeroLag.BinMs;

funcsynapses.ZeroLag.CnxnEndTimesVsRefSpk = synEnds;
funcsynapses.ZeroLag.CnxnEndBinsVsRefSpk = funcsynapses.ZeroLag.CnxnEndTimesVsRefSpk./funcsynapses.ZeroLag.BinMs;
funcsynapses.ZeroLag.CnxnEndTimesVsCCGStart = synEnds+funcsynapses.ZeroLag.BinMs*floor(length(tR)/2);
funcsynapses.ZeroLag.CnxnEndBinsVsCCGStart = ceil(length(tR)/2)+funcsynapses.ZeroLag.CnxnEndTimesVsRefSpk./funcsynapses.ZeroLag.BinMs;

funcsynapses.ZeroLag.CnxnWeightsZ = synWeightsZ;
funcsynapses.ZeroLag.CnxnWeightsR = synWeightsR;
funcsynapses.ZeroLag.PairUpperThreshs = triu(synSig)';
funcsynapses.ZeroLag.PairLowerThreshs = tril(synSig);
% funcsynapses.ZeroLag.TransferWeights = synTransfers;

funcsynapses.ZeroLag.fullCCGMtx = ccgR;%really not very big in terms of bytes
funcsynapses.ZeroLag.CCGbins = tR;
% funcsynapses.CellRates = Rate(S);


