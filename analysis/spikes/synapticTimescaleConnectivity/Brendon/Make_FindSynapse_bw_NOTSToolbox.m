function [funcsynapses] = Make_FindSynapse_bw_NOTSToolbox(basename)

%% Parameters
alpha = 0.01;%signifiance cutoff/pvalue
BinSize = 0.5;%in ms
synapticwindow = [1 4];%Milliseconds in which synapses much occur (inclusive)...
                        % = synWin input to FindSynapses.m
convolutionwidth = 12; % size of the convolution window (in number of bins
                        % = convWin input to FindSynapses.m
excSigWin = 1; % duration of actual ccg above threshold for exc. synapses
inhSigWin = 2; % duration of actual ccg below threshold for inh. synapses
% 
% [T,G] = oneSeries(S);
% T = TimePoints(T);

% SampleRate = 10000;%%Based on usual TSD object sample rate<< NEED TO GENERALIZE THIS, HAVEN'T FIGURED OUT HOW YET

%%%%% added in compensation for no tst
[spiket, spikeind, numclus, iEleClu, ~] = ReadEl4CCG2(basename);
shank = iEleClu(:,2);
cellIx = iEleClu(:,3);
T = spiket;
G = spikeind;
S = 1:numclus;
SampleRate = 20000;
%%%%%

numBinsBinSize = BinSize*SampleRate/1000;
HalfBins = round(300/numBinsBinSize);


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

       if shank(x)~=shank(y)
%            [h,b] = CrossCorr(rgx,rgy,bin,nbBins);
%            cch = [cch h*rx*bin/1000];[funcsynapses] = Make_FindSynapse_bw(S,shank)
           cch(:,end+1) = ccgR(:,x,y);
           cellPairId = [cellPairId;[x y]];
       else
%            [h,b] = CrossCorr(rgx,rgy,bin,nbBinsLg);
%            cchSame = [cchSame h*rx*bin/1000];
%%%%% added in compensation for no tst
           theseindsx = find(spikeind==x);
           theseindsy = find(spikeind==y);
           T = cat(1,spiket(theseindsx),spiket(theseindsy));
           G = cat(1,ones(size(theseindsx)),2*ones(size(theseindsy)));
%%%%%

%            T = cat(1,TimePoints(S{x}),TimePoints(S{y}));
%            G = cat(1,ones(length(S{x}),1),2*ones(length(S{y}),1));
            
           [ccgR_temp, tR_temp, Pairs_temp] = CCG(T, G, numBinsBinSize, HalfBins*10, SampleRate, unique(G), 'count'); %calc cross correlograms, output as counts... 3D output array
           cchSame(:,end+1) = ccgR_temp(:,1,2);
           cellPairIdSame = [cellPairIdSame;[x y]];
       end       
    end
end 
close(h)

[synLat,synStrZ,synStrR,bounds,synStart,synEnd,synFlip] = FindSynapse_bw(cch,'bins',BinSize,'alpha',alpha,'synwin',synapticwindow,'convWin',convolutionwidth);
synTimes = pair2mat(synLat,cellPairId,length(S));
synWeightsZ = pair2mat(synStrZ,cellPairId,length(S));
synWeightsR = pair2mat(synStrR,cellPairId,length(S));
synSig = pair2mat(bounds,cellPairId,length(S));
synStarts = pair2mat(synStart,cellPairId,length(S));
synEnds = pair2mat(synEnd,cellPairId,length(S));
synRev = pair2mat(synFlip,cellPairId,length(S));

[synLat,synStrZ,synStrR,bounds,synStart,synEnd,synFlip] = FindSynapse_SameShank_bw(cchSame,'bins',BinSize,'alpha',0.01,'synwin',synapticwindow,'convWin',convolutionwidth);
synTimes = pair2mat(synLat,cellPairIdSame,synTimes);
synWeightsZ = pair2mat(synStrZ,cellPairIdSame,synWeightsZ);
synWeightsR = pair2mat(synStrR,cellPairIdSame,synWeightsR);
synSig = pair2mat(bounds,cellPairIdSame,synSig);
synStarts = pair2mat(synStart,cellPairIdSame,synStarts);
synEnds = pair2mat(synEnd,cellPairIdSame,synEnds);
synRev = pair2mat(synFlip,cellPairIdSame,synRev);
synRev(isnan(synRev)) = 0;

%% Flipping outputs of above so pre-synaptic cells are on dim 1, Post-syn on dim 2
synTimes = synTimes';
synWeightsZ = synWeightsZ';
synWeightsR = synWeightsR';
synSig = synSig';
synStarts = synStarts';
synEnds = synEnds';
synRev = synRev';

% %% get total connection strengths based on start and end times of each interaction
% [pre,post] = find(~isnan(synWeightsZ));
% strengths = spikeTransfer_InPairSeries(S,[pre post],BinSize,synStarts,synEnds)
% synTransfers = nan(size(synStarts));
% for a = 1:length(pre)
%     synTransfers(pre(a),post(a)) = strengths(a);
% end

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
% funcsynapses.CellRates = Rate(S);
%%%%%%%%%%%%%%%%%%%
totaltime = max(spiket)/SampleRate;
for a = 1:numclus
    id = iEleClu(a,1);
    cellcounts(a) = sum(spikeind==a);
end 
funcsynapses.CellRates = cellcounts/totaltime;
%%%%%%%%%%%%%%%%%%%

