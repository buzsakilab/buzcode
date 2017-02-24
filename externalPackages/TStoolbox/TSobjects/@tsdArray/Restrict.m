function S1 = Restrict(S0,is)

%  Returns a tsdArray restricted in a given intervalSet
%  
%  USAGE:
%  tsO = Restrict(tsa,is) 
%  
%  INPUTS:
%  tsa - a tsdArray
%  is - an intervalSet

%  copyright (c) 2009 Adrien Peyrache adrien.peyrache@gmail.com

S1={};

for i=1:length(S0)
	
	Sr = S0.C{i};
	S1 = [S1;{Restrict(Sr,is)}];

end

S1 = tsdArray(S1);
	