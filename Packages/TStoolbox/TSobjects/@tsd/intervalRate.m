function ic = intervalRate(tsa, is, TimeUnits, varargin)

%  Returns rates of a TSD in each subset of an interval
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
%  		(see OPTIONS for possibilities) and the data gives the rate of
%  		the tsd in each one of the intervals.
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
  
  
  if nargin == 2
    TimeUnits = time_units('s');
  else 
    TimeUnits = time_units(TimeUnits);
  end
    
  opt_varargin = varargin;
  
  time_opt_values = dictArray({ { 'start', []},
		                { 'end', []}, 
		                { 'middle', []} });
  defined_options = dictArray({ { 'time', {'start', {'char'} } } } );
  
  getOpt;

  
  t_ic = intervalCount(tsa, is, 'time', time);


  l = (length(is, TimeUnits));
  if ~isempty(l.data)
    ic = tsd(t_ic.t, t_ic.data ./ l.data);
  else
    ic = 0;
  end
  
  
  