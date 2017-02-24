function ok = belong(is,t)

%  Tells whether timestamps belong to a intervalSet
%  	
%  	USAGE:
%  	ok = belong(is,t) 
%  	
%  	returns 1 if t belongs to the interval is, 0 otherwise
%  	
%  	INPUTS:
%  	is - an interval set
%  	t  - a time expressed in the same Time_Unit, or a vector of times
%  	
%  	OUTPUTS:
%  	ok - a boolean variable, or a vector of booleans

% copyright (c) 2004 Francesco P. Battaglia & Adrien Peyrache
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

  if (size(t,1)~=1 & size(t,2)~=1 ) | length(size(t))>2
	error('t must be a vector');
  end  
  st = Start(is);
  en = End(is);


  ok=zeros(size(t));
  for i=1:length(t)
	for j=1:length(st)
	   if t(i)>=st(j) & t(i)<=en(j)	
		ok(i)=1;
	   end
	end
  end
