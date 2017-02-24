function s = times(a,b)
% s = times(a,b) overload of the .* operator for two tsd objects
%
% the .* operaetor is defined to return a tsd with same timestamps and
% data that are the difference of the data. Returns error if timestamps of the
% two object don't coincide
  
% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

  if ~(isa(a, 'tsd') & isa(b, 'tsd'))
    error(' .* operator defined for two tsd''s!');
  end
  
  if ~compatible_ts(a, b)
    error('Timestamps mismatch')
  end
  
  
  s = tsd(a.t, a.data.*b.data);