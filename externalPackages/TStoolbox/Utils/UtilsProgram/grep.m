function [g, idx]  = grep(C, re)
% [g, idx]  = grep(C, re) similar to UNIX grep function
%   
% INPUTS: 
% C: a cell array of strings
% re: a matlab regular expression
% OUTPUTS:  
% g: the cell array of the strings matched by the regexp
% idx: the indces of the elements in C which matched
  
  
  
  
  [s, e] = regexp(C, re);
  idx = ~cellfun('isempty', s);
  
  for i = 1:length(C)
    if isempty(s{i})
      g{i} = '';
    else
      g{i} = C{i}(s{i}(1):e{i}(1));
    end
  end
  
      