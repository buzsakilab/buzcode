function R = cat(varargin)

%  Concatenates IntervalSet objects
%  
%  	USAGE
%  	R = cat(O1, O2, ..., On, shift)
%  	
%  	INPUTS:
%  	O1, O2, ..., On - intervalSet objects to be concatenated
%  	shift - (optional) a ts object of times each intervalSet object get shifted of
%  		before merging
%  	
%  	OUTPUT:
%  	R - the merged intervalSet object

  
% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
  
  
if nargin < 2
  error('Call like this: R = cat(O1, O2, ..., On, shift)');
end

for i = 1:nargin-1
  if ~isa(varargin{i}, 'intervalSet')
    error('Call like this: R = cat(O1, O2, ..., On, shift)');
  end
end


if isa(varargin{end}, 'intervalSet')
  n_objs = nargin;
  shift = ts(zeros(1, nargin));
elseif length(varargin{end}) == nargin -1 
  n_objs = nargin-1;
else
  error('shift must have as many elements as there are intervalSet  objs to cat');
end

shift = Range(shift, 'ts');

start = zeros(0,1);
stop = zeros(0,1);
for i = 1:n_objs
  O = varargin{i};
  start = [start; Start(O)] + shift(i);
  stop = [stop; End(O)] + shift(i);

end

R = intervalSet(start, stop);