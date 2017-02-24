function [S,ix] = Start(O,tsa)

%  Returns subsets of an intervalSet to which timestamps belong to
%  
%  	USAGE:
%  	S = Restrict(O, TimeUnits)
%  	
%  	INPUTS: 
%  	O   - an intervalSet object
%  	tsa - a tsd or ts object
%  	TimeUnits - a units object or the abbreviation string
%  	
%  	OUTPUT:
%  	S  - an intervalset made of subsets of O which contained a time from tsa.
%  	ix - indices of the subsets containing events from tsa

% copyright (c) 2004 Francesco P. Battaglia, 2008 Adrien Peyrache
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

st = O.start;
ix = [];

for i=1:length(st)

	r = Range(Restrict(tsa,subset(O,i)));
	ix = [ix;length(r)>0];

end

ix2 = find(ix);
S = subset(O,ix2);
