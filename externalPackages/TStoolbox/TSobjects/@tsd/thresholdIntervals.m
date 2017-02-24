function is = thresholdIntervals(tsa, thr, varargin)

%  Returns intervals in which a TSD is above (below) threshold
%  
%  	USAGE:
%  	is = thresholdIntervals(tsa, thr, OptionName, OptionValue)
%  	
%  	INPUTS:
%  	tsa - a tsd object
%  	thr - a threshold value  
%  	
%  	OUTPUTS:
%  	is - the intervalSet of the times in which tsa is above (below) threshold  
%  	
%  	OPTIONS:
%  	'Direction' - possible values are 'Above' (default) 'Below'

% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

 opt_varargin = varargin;
 
 defined_options = dictArray({ { 'Direction', {'Above', { 'char' } } } ...
		   });
 
 getOpt;
 
 switch Direction
  case 'Above'
   st = Range(threshold(tsa, thr, 'Crossing', 'Rising', 'InitialPoint', 1));
   en = Range(threshold(tsa, thr, 'Crossing', 'Falling', 'FinalPoint', 1));
  case 'Below'
   st = Range(threshold(tsa, thr, 'Crossing', 'Falling', 'InitialPoint', 1));
   en = Range(threshold(tsa, thr, 'Crossing', 'Rising', 'FinalPoint', 1));
  otherwise
   error('Unrecognized option value');
 end
 
 
 is = intervalSet(st, en, '-fixit');
 
 
   