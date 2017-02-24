function raster = Raster(S,center,TStart,TEnd, varargin)


opt_varargin = varargin;

defined_options  = dictArray({  { 'BinSize', {10, {'numeric'}}},
                                });
getOpt;

is = intervalSet(Range(center)+TStart, Range(center)+TEnd);
sweeps = intervalSplit(S, is, 'OffsetStart', TStart);

raster = zeros((TEnd-TStart)/BinSize+1,length(sweeps));

for i=1:length(sweeps)

	t = intervalRate2(sweeps{i}, regular_interval(TStart, TEnd, BinSize));
	raster(:,i) = Data(t);
	
end

raster = tsd(Range(t),raster);