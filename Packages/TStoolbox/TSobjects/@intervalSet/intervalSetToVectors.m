function [s, e] = intervalSetToVectors(intervals)

%  De-constructor for intervalSet objects, converts out to two
%  	
%  	USAGE:
%  	[intvect1,intvect2] = intervalSetToVectors(intervals);
%  	
%  	INPUTS:
%  	intervals = any intervalSet object
%
%   OUTPUTS
%   s = vector of start times for collection of intervals
%   e = vector of end times for collection of intervals
% copyright (c) 2014 Brendon O. Watson
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
  
s = Start(intervals);
e = End(intervals);

