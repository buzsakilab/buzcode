function maxTsd = max(tsa,is)

%  Returns the maximal value(s) and time(s) of a TSD object (overoads the max function)
%  
%  	USAGE:
%  	maxTsd = max(tsa,epoch) 
%  	
%  	INPUTS:
%  	tsa - a tsd object
%  	epoch (optional) - an intervalSet object where max values have to be looked for.
%  	
%  	OUTPUTS:
%  	maxTsd - tsd of maximal value(s)
%  
%  copyright (c) 2009 Adrien Peyrache, adrien.peyrache@gmail.com
%  This software is released under the GNU GPL
%  www.gnu.org/copyleft/gpl.html

if nargin==1
	is = timeSpan(tsa);
end

st = Start(is);
en = End(is);

rg = tsa.t;
d = tsa.data;

if size(d,2)>1
	error('The max function for tsd still doesn''t manage multidimensionnal data')
end

ld = length(d);
ls = length(st);
ix=1;

maxVal = [];
maxTimes = [];

for i=1:ls

	while rg(ix)<st(i)
		ix=ix+1;
	end

	ixS = ix;

	while rg(ix)<en(i)
		ix=ix+1;
	end
	
	[m,t] = max(d(ixS:ix));
	maxVal = [maxVal;m];
	maxTimes = [maxTimes;rg(t+ixS-1)];
	
end

maxTsd = tsd(maxTimes,maxVal);