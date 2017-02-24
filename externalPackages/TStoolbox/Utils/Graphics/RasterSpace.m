function spacePETH = RasterSpace(S, phiG, TStart, TEnd, varargin)

opt_varargin = varargin;

defined_options  = dictArray({ { 'RasterFraction', {0.7, {'numeric'}} }
                               { 'BinSize', {10, {'numeric'}}},
                                {'LineWidth', {1, {'numeric'} } },
                                {'Markers', { {}, {'cell'}} } ,
                                 {'MarkerTypes', { {}, {'cell'}}}, 
                                {'MarkerSize', { [], {'numeric'} } }
                                {'SpaceBin', {0.05, {'numeric'} } }

                                });
getOpt;

is = intervalSet(TStart,TEnd);
S = Restrict(S,is);

binnedSpace = [0:SpaceBin:1];
phiR = Restrict(phiG,is);
dPhi = Data(phiR);
spacePETH = zeros(length(binnedSpace),1);

occ = hist(dPhi,binnedSpace);

p = Restrict(phiR,S);
nb = hist(data(p),binnedSpace);
warning off;
spacePETH = nb./occ;
warning on;
spacePETH(isnan(spacePETH))=0;

%  for i=1:length(binnedSpace)-1
%  
%  	ep1 = thresholdIntervals(phiR,binnedSpace(i),'Direction','Above');
%  	ep2 = thresholdIntervals(phiR,binnedSpace(i+1),'Direction','Below');
%  	ep = intersect(ep1,ep2);
%  	spacePETH(i) = nanmean(Data(intervalRate(S,ep)));
%  
%  end

