function O = realign(tsa, varargin)

%  Realign timestamps  
%  	
%  	USAGE
%  	tso = realign(tsa, OptionName, OptionValue) 
%  	
%  	INPUTS:
%  	O: a tsd object
%  	
%  	OUTPUTS: 
%  	O: the realigned tsd  
%  	
%  	OPTIONS:
%  	'ZeroFirst' - if set to non-zero, the timestamp will be realigned so
%  		that the first timestamp corresponds to zero  
%  	'ZeroLast'  - if set to non-zero, the timestamp will be realigned so
%  		that the last timestamp corresponds to zero  
%  	'ZeroAt'    - if set to non-zero, the timestamp will be realigned so
%  		that the specified value will correspond to zero 

% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
  
  
  opt_varargin = varargin;

  defined_options = dictArray({ { 'ZeroFirst', {0, {'numeric'} } }, 
		                { 'ZeroLast', {0, {'numeric'} } },
		                { 'ZeroAt', {0, {'numeric'} } } });
  
  getOpt;
  
  if ((ZeroFirst ~= 0) + (ZeroLast ~= 0) + (ZeroAt ~= 0)) ~= 1
    error('Set exactly one option among ZeroFirst, ZeroLast, ZeroAt');
  end
  
  
  t = tsa.t;
  
  if ZeroFirst ~= 0
    t = t - t(1);
  elseif ZeroLast ~= 0
    t = t - t(end);
  elseif ZeroAt
    t = t - ZeroAt;
  end
  
  O = tsd(t, tsa.data);
  
    