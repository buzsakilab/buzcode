function EnO = End(S0)

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

EnO = [];

for i=1:length(S0)
	
	Sr = S0.C{i};
	rg = Range(Sr);
	if length(rg)>0
		EnO = [EnO;rg(end)];
	end

end
