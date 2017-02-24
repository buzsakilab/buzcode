function S = smooth(tsa, l, varargin)

%  Smooth a TSD
%  	
%  	USAGE:
%  	S = smooth(tsa, l, OptionName, OptionValue, ...) 
%  	%
%  	This is intended for continuously sampled tsd-s, the sampling rate is
%  	guessed with the median inter-event interval method (see MEDIAN_DT)
%  	If there are some missing values, they are interpolated before smoothing
%  	
%  	INPUTS:
%  	tsa - a tsd object
%  	l   - the length of the smoothing window, in the units specified by the
%  	TimeUnits option
%  	
%  	OUTPUTS:
%  	S - the smoothed tsd 
%  	
%  	OPTIONS:
%  	'TimeUnits' - specifies the time units for l (defaults to tsa.time_unit)
%  	'UseWindow' - type of window to use for smoothing, defaults to
%  		hamming. Admissible   values are 'bartlett', 'blackman', 'boxcar',
%  		'chebwin', 'hamming', 'hann', 'kaiser'
%  		'InIntervals', if non-empty, must be a cell array of intervalSet.
%  		In that case, the tsd is taken to be non-contiguous, 
%  		and to only assume values in the specified set of IntervalSets
 
% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
 
  
  opt_varargin = varargin;
  defined_options = dictArray( ...
  { { 'TimeUnits', {tsa.time_unit, {'char', 'units'} } }, 
    { 'UseWindow', {'hamming', {'char'} } },  
    {'InIntervals', { {}, {'cell'} } } 
    });
  
  getOpt;
  
  if ~isempty(InIntervals)
      for i = 1:length(InIntervals)
          2+3;
          
          tsdR{i} = Restrict(tsa, InIntervals{i});
          tsdR{i} = smooth(tsdR{i}, l, 'TimeUnits', TimeUnits, 'UseWindow', UseWindow);
      end
      S = cat(tsdR{:});
      return
  end
  
  l = l / median_dt(tsa, TimeUnits);
  
  l = floor(l/2) * 2;
  
  eval(['hh = ' UseWindow '(l, ''symmetric'');']);
  
  tsi = interp_missing(tsa);
  
  v = Data(tsi);

  t = Range(tsi);
  
  v = conv(v, hh);
  
  v = v(l/2:end-l/2);

  v = v / sum(hh);

  S = tsd(t, v);
