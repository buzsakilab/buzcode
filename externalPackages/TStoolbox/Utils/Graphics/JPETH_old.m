function [M,dt] = JPETH(T1,T2,Tref,bin,win,varargin)

% USAGE:
% 	[C,B] = JPETH(T1,T2,Tref,bin,win,options)
% 
% T1,T2 : tsd from whom JPETH are drawn
% TRef : tsd represneting the common reference for T1,T2
% bin : larger of bins (in ms)
% win : numbers of bins
% Option:
%	'LastSpk' : take for T2 only the nearest events from TRef. 
%
% Adrien Peyrache 2007

%  keyboard

if length(varargin)==1
	if varargin{1} == 'LastSpk';
		ls = 1;
	else
		error(['Option ' varargin{1} ' not suppoerted']);
	end

elseif length(varargin)>1
	error(['Too many options']);
else
	ls = 0;
end

M = zeros(win+1,win+1);
bin = bin*10;

warning off
dt = [-win/2*bin:bin:win/2*bin];

is = intervalSet(Range(Tref)-win/2*bin, Range(Tref)+win/2*bin);
sweeps1 = intervalSplit(T1, is, 'OffsetStart', -win/2*bin);
sweeps2 = intervalSplit(T2, is, 'OffsetStart', -win/2*bin);

for i=1:length(Start(is))

	C1=hist(Range(sweeps1{i}),dt);
	C2=hist(Range(sweeps2{i}),dt);

	if ls & sum(C2)>0

		ix = find(C2);
		[d , pos] = min(abs(ix - win/2));
		C2t = zeros(size(C2));
		C2t(ix(pos)) = C2(ix(pos));
%  		keyboard
		C2 = C2t;
	
	end

	t = C2'*C1;
	M = M + t;

end

M = 10000*M/(length(Range(Tref))*bin);
dt = dt/10;
