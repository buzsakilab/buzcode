function is = timeSpan(S)

%  Returns interval spanning the duration of the tsdArray
%  	
%  	USAGE:
%  	is = timeSpan(tsa)
%  	
%  	OUTPUT:
%  	is - an intervalSet object
%  
%  copyright (c) 2009 Adrien Peyrache adrien.peyrache@gmail.com

T_start = inf;
T_end = -inf;

for iC = 1:length(S.C)
	if ~isempty(Data(S.C{iC}))
	T_start = min(T_start, StartTime(S.C{iC}));
	T_end = max(T_end, EndTime(S.C{iC}));
	end
end

is = intervalSet(T_start,T_end);