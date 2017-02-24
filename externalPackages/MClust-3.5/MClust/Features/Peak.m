function peak = Peak(WV)
% p2v = Peak(WV)%% INPUTS:%     WV = tsd of nSpikes x nTrodes (4) x nChannels%% OUTPUTS:%       p2v = nSpikes x nTrodes (4) x spike width%% ALGO:%       ratio of peak height:valley depth%% ADR 1998%% Status: PROMOTED (Release version) % See documentation for copyright (owned by original authors) and warranties (none!).% This code released as part of MClust 3.0.% Version control M3.0.
WVD = Data(WV);peak = max(WVD, [], 3);S = warning;warning offwarning(S);
