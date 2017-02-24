function ic = intervalCount(tsa, is, varargin)

%  Counts instances of tsa in each subset of an interval
%  	
%  	USAGE:
%  	ic = intervalCount(tsa, is, options) 
%  	
%  	INPUTS:
%  	tsa 	- a tsd object
%  	is	- an intervalSet
%  	
%  	OUTPUTS: 
%  	ic 	- a tsd object, where the timestamps correspond to each interval 
%  		(see OPTIONS for possibilities) and the data gives the number of points
%  		in the tsd in each one of the intervals
%  	
%  	OPTIONS:
%  	'time' determines which time is selected for each interval, possible
%  		values are :
%  		'start'  - use start of intervals (default)
%  		'end'    - use end of intervals
%  		'middle' - use middle point of intervals

% copyright (c) 2004 Francesco P. Battaglia, 2009 Adrien Peyrache adrien.peyrache@gmail.com
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

  opt_varargin = varargin;
  
  
  time_opt_values = dictArray({ { 'start', []},
		                { 'end', []}, 
		                { 'middle', []} });
  defined_options = dictArray({ { 'time', {'start', {'char'} } } } );
  
  getOpt;

st = Start(is, tsa.time_unit);
en = End(is,tsa.time_unit);

spk = tsa.t;
l = length(spk);
ic = zeros(length(st),1);

ix = 1; %index of spike times. Begins at 1!
if length(spk)
for j=1:length(st) % runs through all bins
	
	% the idea in this piece of code is to scan all spikes
	% and check if it is in the given bin.

	while spk(ix)<st(j) & ix<l 
		ix = ix+1;
	end
	while spk(ix)<=en(j) & ix<l 
		ic(j) = ic(j)+1;
		ix = ix+1;
	end
	if ix==l & spk(ix)<=en(j) & spk(ix)>st(j) 
		ic(j) = ic(j)+1;
	end

end
end

switch time
	case 'start'
	t_ic = Start(is, tsa.time_unit);
	case 'end'
	t_ic = End(is, tsa.time_unit);
	case 'middle'
	t_ic = ( Start(is, tsa.time_unit) + End(is, tsa.time_unit) ) / 2;
end
  
ic = tsd(t_ic, ic);