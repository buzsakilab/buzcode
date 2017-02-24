function is = boundaryIntervals(tsa, thrLow,thrHigh)

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
 if any(Data(tsa)<thrLow) 
 	epLow = thresholdIntervals(tsa,thrLow,'Direction','Above');
 else
 	epLow = timeSpan(tsa);
 end
 if any(Data(tsa)>thrHigh) 
 	epHigh = thresholdIntervals(tsa,thrHigh,'Direction','Below');
 else
 	epHigh = timeSpan(tsa);
 end

 is = intersect(epLow,epHigh);
 
 
   