function ix = findIdx(C, TSA);
%  ix = findIdx(C, TSA) find  compound events in a tsd, return array of
%  indices of event centers
% INPUTS:
% C: a compoundEvents object
% TSA: a tsd object
% OUTPUTS:
% ix: the array of the indices into the TSA of the event centers

% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
  
  % keyboard
  ix = findIdx_c(C, Range(TSA));
  ix = ix + C.condOffset;
  
  