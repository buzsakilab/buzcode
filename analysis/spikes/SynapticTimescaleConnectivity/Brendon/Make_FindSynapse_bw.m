function [funcsynapses] = Make_FindSynapse_bw(S,shank,cellIx)

%% Parameters
alpha = 0.01;%signifiance cutoff/pvalue
BinSize = 0.5;%in ms
synapticwindow = [1 4];%Milliseconds in which synapses much occur (inclusive)...
                        % = synWin input to FindSynapses.m
convolutionwidth = 12; % size of the convolution window (in number of bins
                        % = convWin input to FindSynapses.m
excSigWin = 1; % duration of actual ccg above threshold for exc. synapses
inhSigWin = 1.5; % duration of actual ccg below threshold for inh. synapses

[T,G] = oneSeries(S);
T = TimePoints(T);

SampleRate = 10000;%%Based on usual TSD object sample rate<< NEED TO GENERALIZE THIS, HAVEN'T FIGURED OUT HOW YET
numBinsBinSize = BinSize*SampleRate/1000;
HalfBins = round(300/numBinsBinSize);

% msextractedperspike = 1.6;%based on typical 32 points at 20000khz

%% Get raw CCGs for all pairs (will not use these for same-shank cells
[ccgR, tR, Pairs] = CCG(T, G, numBinsBinSize, HalfBins, SampleRate, unique(G), 'count'); %calc cross correlograms, output as counts... 3D output array

%% Iterate through all pairs
cch = [];
cchSame = [];

cellPairId = [];
cellPairIdSame = [];

for x=1:length(S)
    h=waitbar(x/length(S));
%     rgx = Range(S{x});
%     rx = length(Range(S{x}));
    for y=x+1:length(S)
%        rgy = Range(S{y});

%        if [x == 55 && y == 62] || [x == 54 && y == 55]
%            1;
%        end
       if shank(x)~=shank(y)
%            [h,b] = CrossCorr(rgx,rgy,bin,nbBins);
%            cch = [cch h*rx*bin/1000];[funcsynapses] = Make_FindSynapse_bw(S,shank)
           cch(:,end+1) = ccgR(:,x,y);
           cellPairId = [cellPairId;[x y]];
       else
%            [h,b] = CrossCorr(rgx,rgy,bin,nbBinsLg);
%            cchSame = [cchSame h*rx*bin/1000];
           T = cat(1,TimePoints(S{x}),TimePoints(S{y}));
           G = cat(1,ones(length(S{x}),1),2*ones(length(S{y}),1));
            
           [ccgR_temp, tR_temp, Pairs_temp] = CCG(T, G, numBinsBinSize, HalfBins*10, SampleRate, unique(G), 'count'); %calc cross correlograms, output as counts... 3D output array
           cchSame(:,end+1) = ccgR_temp(:,1,2);
           cellPairIdSame = [cellPairIdSame;[x y]];
       end       
    end
end 
close(h)

[synLat,synStrZ,synStrR,bounds,synStart,synEnd,synFlip] = FindSynapse_bw(cch,'bins',BinSize,'alpha',alpha,'synwin',synapticwindow,'convWin',convolutionwidth,'excSigWin',excSigWin,'inhSigWin',inhSigWin);
synTimes = pair2mat(synLat,cellPairId,length(S));
synWeightsZ = pair2mat(synStrZ,cellPairId,length(S));
synWeightsR = pair2mat(synStrR,cellPairId,length(S));
synSig = pair2mat(bounds,cellPairId,length(S));
synStarts = pair2mat(synStart,cellPairId,length(S));
synEnds = pair2mat(synEnd,cellPairId,length(S));
synRev = pair2mat(synFlip,cellPairId,length(S));

[synLat,synStrZ,synStrR,bounds,synStart,synEnd,synFlip] = FindSynapse_SameShank_bw(cchSame,'bins',BinSize,'alpha',0.01,'synwin',synapticwindow,'convWin',convolutionwidth,'excSigWin',excSigWin,'inhSigWin',inhSigWin);
synTimes = pair2mat(synLat,cellPairIdSame,synTimes);
synWeightsZ = pair2mat(synStrZ,cellPairIdSame,synWeightsZ);
synWeightsR = pair2mat(synStrR,cellPairIdSame,synWeightsR);
synSig = pair2mat(bounds,cellPairIdSame,synSig);
synStarts = pair2mat(synStart,cellPairIdSame,synStarts);
synEnds = pair2mat(synEnd,cellPairIdSame,synEnds);
synRev = pair2mat(synFlip,cellPairIdSame,synRev);
synRev(isnan(synRev)) = 0;

%%

%% Flipping outputs of above so pre-synaptic cells are on dim 1, Post-syn on dim 2
synTimes = synTimes';
synWeightsZ = synWeightsZ';
synWeightsR = synWeightsR';
synSig = synSig';
synStarts = synStarts';
synEnds = synEnds';
synRev = synRev';

%% get total connection strengths based on start and end times of each interaction
% [pre,post] = find(~isnan(synWeightsZ));
% strengths = spikeTransfer_InPairSeries(S,[pre post],funcsynapses,synStarts,synEnds);
% synTransfers = nan(size(synStarts));
% for a = 1:length(pre)
%     synTransfers(pre(a),post(a)) = strengths(a);
% end

% %% Linearly scale any bins partially cut by the width of spike sampling
% % will not be fully accurate since spiking dynamics will not be equal
% % across the bin
% binstarts = tR-BinSize/2;%tR gives bin centers, this gives actual start times...
% binends = tR+BinSize/2;%...and stop times
% cutbins(1) = find(binstarts<-msextractedperspike/2 & binends>-msextractedperspike/2);%figure out which bins have been cut
% cutbins(2) = find(binstarts<msextractedperspike/2 & binends>msextractedperspike/2);
% remprop = rem(msextractedperspike,floor(msextractedperspike))/BinSize;%proportion of usual bin width occupied by the cut bins
% convertor = 1/remprop;
% ScaledBins = ccgR(cutbins,:,:)*convertor;

%% Prep for output
funcsynapses.BinMs = BinSize;
funcsynapses.AlphaValue = alpha;
funcsynapses.SynapticWindowInMs = synapticwindow;
funcsynapses.ConvolutionWidthInBins = convolutionwidth;
funcsynapses.ExcMinimumWidth = excSigWin;
funcsynapses.InhMinimumWidth = inhSigWin;
funcsynapses.CellShanks = shank;
funcsynapses.CellShankIDs = cellIx;
funcsynapses.FlippedCnxns = synRev;

funcsynapses.CnxnTimesVsRefSpk = synTimes;
funcsynapses.CnxnBinsVsRefSpk = funcsynapses.CnxnTimesVsRefSpk./funcsynapses.BinMs;
funcsynapses.CnxnTimesVsCCGStart = synTimes+funcsynapses.BinMs*floor(length(tR)/2);
funcsynapses.CnxnBinsVsCCGStart = ceil(length(tR)/2)+funcsynapses.CnxnTimesVsRefSpk./funcsynapses.BinMs;

funcsynapses.CnxnStartTimesVsRefSpk = synStarts;
funcsynapses.CnxnStartBinsVsRefSpk = funcsynapses.CnxnStartTimesVsRefSpk./funcsynapses.BinMs;
funcsynapses.CnxnStartTimesVsCCGStart = synStarts+funcsynapses.BinMs*floor(length(tR)/2);
funcsynapses.CnxnStartBinsVsCCGStart = ceil(length(tR)/2)+funcsynapses.CnxnStartTimesVsRefSpk./funcsynapses.BinMs;

funcsynapses.CnxnEndTimesVsRefSpk = synEnds;
funcsynapses.CnxnEndBinsVsRefSpk = funcsynapses.CnxnEndTimesVsRefSpk./funcsynapses.BinMs;
funcsynapses.CnxnEndTimesVsCCGStart = synEnds+funcsynapses.BinMs*floor(length(tR)/2);
funcsynapses.CnxnEndBinsVsCCGStart = ceil(length(tR)/2)+funcsynapses.CnxnEndTimesVsRefSpk./funcsynapses.BinMs;

funcsynapses.CnxnWeightsZ = synWeightsZ;
funcsynapses.CnxnWeightsR = synWeightsR;
funcsynapses.PairUpperThreshs = triu(synSig)';
funcsynapses.PairLowerThreshs = tril(synSig);
% funcsynapses.TransferWeights = synTransfers;

funcsynapses.fullCCGMtx = ccgR;
funcsynapses.CCGbins = tR;
% funcsynapses.CutBinIdxs = cutbins;
% funcsynapses.TotalCutWidthAroundSpike = msextractedperspike;
% funcsynapses.CutBinConversion = convertor;
% funcsynapses.ConvertedBins = ScaledBins;
funcsynapses.CellRates = Rate(S);


