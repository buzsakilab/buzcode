function [ crange ] = SpecColorRange( spec,range )
%[ crange ] = SpecColorRange( spec,range ) calculates a good color range for a
%spectrogram and sets the plot to the good range.
%
%
%DLevenstein 2017
%%
if ~exist('range','var')
    range = [2.5 2.5];
end

[~,mu,sig] = zscore(spec);
crange = [min(mu)-range(1)*max(sig) max(mu)+range(2)*max(sig)];
caxis(crange)

end

