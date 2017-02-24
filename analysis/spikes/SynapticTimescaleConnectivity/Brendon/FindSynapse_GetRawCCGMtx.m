function funcsynapses = FindSynapse_GetRawCCGMtx(funcsynapses,rawS)
% 

%% Parameters
SampleRate = 10000;
BinSize = funcsynapses.BinMs;
numBinsBinSize = BinSize*SampleRate/1000;
HalfBins = round(300/numBinsBinSize);

%% Convert spiking info format 
[T,G] = oneSeries(rawS);
T = TimePoints(T);

%% Get raw CCGs for all pairs (will not use these for same-shank cells
[ccgR, dummy, dummy] = CCG(T, G, numBinsBinSize, HalfBins, SampleRate, unique(G), 'count'); %calc cross correlograms, output as counts... 3D output array

%% save raw ccg
funcsynapses.fullRawCCGMtx = ccgR;
