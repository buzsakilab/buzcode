function T1 = StartTime(tsa, tsflag)

% Returns first timestamp
%  	
%  	USAGE:
%  	T1 = StartTime(tsa, tsflag)
%  	
%  	OPTIONS:
%  	tsflag - if 'ts' returns time in timestamps (default),
%  		 if 'sec' returns time in sec
%  		 if 'ms' returns time in ms

% the current version is an evolution of A. David Redish's code
% changes copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% released under the GPL   

  
T1 = tsa.t(1);

units_out = tsa.time_unit; 



if nargin == 2
  if isa(tsflag, 'char')
    units_out = time_units(tsflag);
  elseif isa(tsflag, 'units')
    units_out = tsflag;
  else
    error(['tsflag must be units object or corresponding abbreviation', ...
	   ' string']);
  end
end

cnvrt = convert(tsa.time_unit, units_out);
if cnvrt ~= 1
  T1 = cnvrt * T1;
end
