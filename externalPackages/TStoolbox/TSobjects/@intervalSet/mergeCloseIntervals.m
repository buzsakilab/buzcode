function iso = mergeCloseIntervals(is, thr, varargin)

%  Merge together intervals that are closer than a given threshold
%  
%  	USAGE:
%  	iso = mergeCloseIntervals(is, thr, OptionName, OptionValue) 
%  	
%  	INPUTS:
%  	is  - an intervalSet
%  	thr - a distance threshold; if two consecutive intervals are closer than
%  	      thr they will be merged 
%  	
%  	OUTPUTS:
%  	iso - the resulting intervalSet
  
% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

  
  
  opt_varargin = varargin;
  
  defined_options = dictArray( ...
  { { 'TimeUnits', {time_units('ts'), {'char', 'units'} } } } );
  
  getOpt;  
  
  
  cnvrt = convert(TimeUnits, is.units);
  
  thr = thr * cnvrt;

if length(is.start)>0
  
  tsp = timeSpan(is);
  
  df = tsp- is;
  
  df = dropShortIntervals(df, thr);
  
  iso = tsp - df;
  
else

  iso = is;

end