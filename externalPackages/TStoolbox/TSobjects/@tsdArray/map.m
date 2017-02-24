function SO = map(S, cmd_string, varargin)

% Applies function on each element of the array
%  
%  	USAGE:
%  	SO = map(S, cmd_string, arg1, ..., argN)
% 
% Similarly to Python's function, map applies the commands indicated in
% cmd_stirng to each elements of the tsdArray S. In cmd_string, the input
% tsd array must be referred to as TSA, whereas the output tsd must be
% referred to as TSO. Further arguments may be used in cmd_string, and
% may be passed to map as arguments following cmd_string. In cmd_string,
% they must be referred to %1, %2, etc. Admittedly a dirty hack, but the
% only way to mimick passing a function object... 
  
% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html  
  
  expand_cmd_string
  
  C = cell(size(S.C));
  
  for i = 1:length(S)
    TSA = S.C{i};
    eval(cmd_string);
    TSO = setName(TSO, Name(TSA));
    C{i} = TSO;
  end
  
  SO = tsdArray(C);
  
     