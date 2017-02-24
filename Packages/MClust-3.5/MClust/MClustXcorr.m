function [Y,xdim] = MClustXcorr(a, b, binsize_msec, xcorr_width_msec)
%function [Y,xdim] = MClustXcorr(cl1, cl2, binsize_msec, msec_each_side)
% INPUT: cluster number 1
%        cluster number 2
%        bin size of choice
%        msec on each side of the xcorr plot.
% OUTPUT:
%  an xcorr plot
%  the y and x dimensions of the plot.
%
% TO USE WITH MCLUST, put this in the MClust/ClusterOptions folder
% cowen
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

global MClust_TTData MClust_Clusters 

ts1 = Range(ExtractCluster(MClust_TTData, FindInCluster(MClust_Clusters{a})),'ts');
ts2 = Range(ExtractCluster(MClust_TTData, FindInCluster(MClust_Clusters{b})),'ts');
nbins = round(xcorr_width_msec/binsize_msec);
if isempty(ts1) || isempty(ts2)
    msgbox('No spikes in one of the chosen clusters.')
else
    [Y, xdim] = CrossCorr(ts1, ts2, binsize_msec, nbins);
    plot(xdim,Y)
    xlabel('msec')
    ylabel('Count')
    title(['XCorr ' num2str(a) ' Vs ' num2str(b) ]);
end

