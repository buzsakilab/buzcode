function SO = JitterPattern(S,bound,varargin)

%  SO = JitterPattern(ts,bound,epoch) returns a ts object SO whose timestamps in the intervalSet epoch are shifted by a jitter belonging to the interval [-bound bound];
%  INPUT:
%  S : a ts object or a tsdArray
%  bound : absolute bound for the jitters in the timestamps unit
%  optionnal:
%  epoch : the intervalSet in which the jittering has to take place.
%  
%  OUTPUT :
%  a ts object or a tsdArray 
%  
%  
%  Adrien Peyrache 2007


if length(varargin) == 1
	if isa(varargin{1},'intervalSet')
		epoch = varargin{1};
		S = Restrict(S,epoch);
	else
		warning('argument must be an intervalSet, not taken into account')
	end
elseif length(varargin)>1
	warning('too many arguments!!')
else
	epoch = intervalSet(min(Start(S)),max(End(S)));

end


if isa(S,'tsdArray')

	SO = {};

	for i=1:length(S)
		SO = [SO;{shift(S{i},bound)}];
	end

	SO = Restrict(tsdArray(SO),epoch);

else

	SO = Restrict(shift(S,bound),epoch);

end

end

function SO = shift(S,bound)

	rg = Range(S);
	jitter = 2*bound*rand - bound;
	rg = rg - jitter;
	SO = ts(rg);

end
	


