function [synTimes,synWeightsZ,synWeightsR,synSig] = Make_FindSynapse(S,shank)

%Parameters
bin = 0.5;
nbBins = 100;
nbBinsLg = 1000;

cch = [];
cchSame = [];

cellPairId = [];
cellPairIdSame = [];

for x=1:length(S)
    h=waitbar(x/length(S));
    rgx = Range(S{x});
    rx = length(Range(S{x}));
    for y=x+1:length(S)
        
       rgy = Range(S{y});
       if shank(x)~=shank(y)
           [h,b] = CrossCorr(rgx,rgy,bin,nbBins);
           cch = [cch h*rx*bin/1000];
           cellPairId = [cellPairId;[x y]];
       else
           [h,b] = CrossCorr(rgx,rgy,bin,nbBinsLg);
           cchSame = [cchSame h*rx*bin/1000];
           cellPairIdSame = [cellPairIdSame;[x y]];
       end       
    end
end
close(h)

[synLat,synStrZ,synStrR,bounds] = FindSynapse(cch,'synWin',[0 8]);
synTimes = pair2mat(synLat,cellPairId,length(S));
synWeightsZ = pair2mat(synStrZ,cellPairId,length(S));
synWeightsR = pair2mat(synStrR,cellPairId,length(S));
synSig = pair2mat(bounds,cellPairId,length(S));

[synLat,synStrZ,synStrR,bounds] = FindSynapse_SameShank(cchSame,'synWin',[0 8]);
synTimes = pair2mat(synLat,cellPairIdSame,synTimes);
synWeightsZ = pair2mat(synStrZ,cellPairIdSame,synWeightsZ);
synWeightsR = pair2mat(synStrR,cellPairIdSame,synWeightsR);
synSig = pair2mat(bounds,cellPairIdSame,synSig);
