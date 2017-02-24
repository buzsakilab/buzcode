function ControlIntervals = ControlIntervals(int1,int2,n)
%  
%  USAGE : 
%  	ControlIntervals = ControlIntervals(int1,int2,n)
%  	this function split an interval int1 in n subintervals whose length is equal to the total length of interval int2
%  
%  Input :
%  	int1 : interval to split
%  	int2 : interval to avoid
%  	n : number of subintervals to compute
%  
%  Adrien Peyrache 2007


T = sum(End(int2)-Start(int2));

tb = Start(int1);
te = End(int1);
l = te(end)-tb(1);

int1 = int1-int2;
st = Start(int1);
en = End(int1);

ControlIntervals = cell(n,1);

for i=1:n

	ok=0;
	while ~ok
	
		t0 = l*rand+tb(1);
		while ~belong(int1,t0)
			t0 = l*rand+tb(1);
		end
	
		ix = 1;
		stCtl = [];
		enCtl = [];
	
		while st(ix)<t0 & ix<length(st)
			ix=ix+1;
		end

		ix = ix-1;
		stCtl(1) = t0;
		enCtl(1) = en(ix); 

		ixCtl = 1;
	
		while sum(enCtl-stCtl)<T & (ix+ixCtl)<length(st)
			stCtl(ixCtl+1) = st(ix+ixCtl);
			enCtl(ixCtl+1) = en(ix+ixCtl);
			ixCtl = ixCtl+1;
		end

		dt = sum(enCtl-stCtl)-T;
		if dt>0
			enCtl(end) = stCtl(end)+dt;
			if enCtl(end)<en(end)
				ok=1;
			end
		end

	end

	ControlIntervals{i} = intervalSet(stCtl,enCtl);

end