function t = threshold(tsa, thr, varargin)

%  Threshold crossings
%  	
%  	USAGE:
%  	t = threshold(tsa, OptionName, OptionValue)
%  	
%  	INPUTS:
%  	tsa - a tsd object
%  	thr - a threshold value 
%  	
%  	OPTIONS:
%  	'Crossing' - type of crossing, possible values are:
%  			'Rising' (default)
%  			'Falling' 

% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

 opt_varargin = varargin;
 
 defined_options = dictArray({ { 'Crossing', { 'Rising' , {'char' } } }, ...
      { 'InitialPoint', {0, {'numeric'} } }  ...
      { 'FinalPoint', {0, {'numeric'} } } }  ...
			     );
 
 getOpt;
 
 d = tsa.data;

 d1 = d(1:(end-1));
 d2 = d(2:end);
 
 
if length(d1)*length(d2)>0
	
	switch Crossing
	case 'Rising'
	
	ix = find(d1 < thr & d2 >= thr) + 1;
    
	if InitialPoint
	if d1(1) > thr
		ix = [1; ix];
	end
    end
    if FinalPoint
	if d1(end) <= thr
		ix = [ix; length(d)];
	end
	end
	case 'Falling'
	ix = find(d1 >= thr & d2 < thr) + 1;
	if InitialPoint
	if d1(1) <= thr
		ix = [1; ix];
	end
    end
    if FinalPoint
	if d1(end) >= thr
		ix = [ix; length(d)];
	end
	end
	otherwise
	error('Unrecognized option value');
	end
	
	
	t = subset(tsa, ix);
 
else
	
	t = tsd([],[]);
end