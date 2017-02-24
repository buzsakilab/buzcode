function r = Rate(S0,varargin)

%  Returns a vector containing the time of the last event of each tsd
%  	
%  	USAGE:
%  	eno = Start(tsa,is)
%  	
%  	INPUTS:
%  	tsa - a tsdArray
%  	is - an intervalSet
%  	
%  	OUTPUT:
%  	eno - an vector of time

%  copyright (c) 2009 Adrien Peyrache adrien.peyrache@gmail.com

r = [];
if nargin>1
    ep = varargin{1};
    S0 = Restrict(S0,ep);    
else
    ep = intervalSet(min(Start(S0)),max(End(S0)));
end
T = tot_length(ep,'s');
r = zeros(length(S0),1);

for c=1:length(S0)
	
	Sr = S0.C{c};
	rg = Range(Sr);
	if ~isempty(rg)
		r(c) = length(rg)/T;
	end

end
