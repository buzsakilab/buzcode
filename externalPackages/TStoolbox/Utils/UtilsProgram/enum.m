function enum(S, start)
% enum(S) fake a C-like enum
%
% enum takes a cell array of strings and assigns variables with their
% names in the caller's workspace. The variables are given increasing integer
% variables starting  1 or from start, if given
% INPUTS: 
% S: a cell aray of strings: variable names 
% start (optional): value of the first variable 
  
% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

  init = 1;
  if nargin > 1
    start = init;
  end
  
  
  for i = 1:length(S)
    assignin('caller', S{i}, init - 1 + i);
  end
  
  