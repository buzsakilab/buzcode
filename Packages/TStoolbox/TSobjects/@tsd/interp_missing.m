function tsi = interp_missing(tsa, TimeUnits, varargin)

%  Interpolate missing values in continuously sampled tsd
%  
%  	USAGE:
%  	tsi = interp_missing(tsa, OptionName, OptionValue, ...) 
%  	
%  	INPUTS:
%  	tsa - a tsd object
%  	
%  	OUTPUTS:
%  	tsi - a tsd object, with the missing values interpolated, assuming the
%  	      sampling rate as the inverse of the median inter-event interval
%  	
%  	OPTIONS:
%  	'InterpNaN' - if set to non-zero, the missing values are set to NaN

% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

  
  
  if nargin == 1
    TimeUnits = time_units('ts');
  else 
    TimeUnits = units(TimeUnits);
  end
  
  opt_varargin = varargin;
  
  defined_options = dictArray({ { 'InterpNaN', {0, {'numeric'} } } } );
  
  getOpt;
  
  [mdt, n_missing] = median_dt(tsa, TimeUnits);
  mdt_min = mean(diff(Range(tsa)));
  if mdt_min < mdt
      warning('strange tail in inter-frame intervals');
      mdt_min = mdt_min / 2;
  else
      mdt_min = mdt * 0.9;
  end
  if n_missing / length(Range(tsa)) > 0.05
    warning('more than 5% of points missing');
  end
  
  t = Range(tsa, TimeUnits);
  
  [nt, orig_ix] = interp_missing_c(t, mdt, mdt_min);
  nd = interp1(t, Data(tsa), nt);

  if InterpNaN
    m_ix = setdiff(1:length(t), orig_ix);
    nd(m_ix) = NaN;
  end
  
  
  tsi = tsd(nt, nd);
  