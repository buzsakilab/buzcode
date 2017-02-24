function is = toIntervalSet(tsa)

%  Returns the intervalSet delimited by the timestamps f a TSD
%  
%  	USAGE:
%  	is = toIntervalSet(tsa) 
%  	
%  	INPUTS:
%  	tsa - a tsd object
%  	
%  	OUTPUTS:
%  	is - the intervalSet where each interval is of the form [tsd(i) tsd(i+1)]


% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

  
t = tsa.t;

is = intervalSet(t(1:end-1), t(2:end));

