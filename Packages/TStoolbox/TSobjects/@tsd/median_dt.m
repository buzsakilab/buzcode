function [mdt, n_missing] = median_dt(tsa, TimeUnits)

%  Find the median value of inter-event interval 
%  	
%  	USAGE:
%  	[mdt [, n_missing] ] = median_dt(tsa) 
%  	
%  	Intended for tsd that are sampled at regular intervals, this function
%  	returns the mean inter-event interval.
%  	
%  	INPUTS:
%  	tsa 	  - a tsd object 
%  	TimeUnits - a units object or the abbreviation string
%  	
%  	OUTPUTS:
%  	mdt	  - the median inter-event interval, in the units specified by TimeUnits 
%  	n_missing - (optional) the number of points that are missing, assuming
%  		    that everything was sampled at 1/mdt rate
  
% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

  
  
  
  if nargin == 1
    TimeUnits = tsa.time_units;
  else 
    TimeUnits = time_units(TimeUnits);
  end
  
  
  dt = diff(Range(tsa, TimeUnits));
  
  mdt = median(dt);
  
  miss_ix = find((dt > mdt + eps ) | (dt < mdt - eps));
  
  if nargout > 1
    n_missing = sum(dt(miss_ix)) / mdt - length(miss_ix);
  end
  
  
    