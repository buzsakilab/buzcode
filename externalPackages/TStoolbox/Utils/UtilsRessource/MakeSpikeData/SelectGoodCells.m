function A = truc(A)

A = getResource(A,'SpikeData');
A = registerResource(A, 'GoodCells', 'numeric', {[], []}, ...
    'gc', ['good cells']);

nbC = length(S);
gc = zeros(nbC,1);


for i=1:nbC

	rg = Range(S{i});
	if length(rg)>0
		[C,B] = CrossCorr(rg,rg,0.1,100);
		t = find((B<3&B>1) | (B>-3 & B<-1));

		T = timeSpan(S{i});
		T = (End(T)-Start(T))/10000;
		asym = length(rg)/T;
		m = sum(C(t)>asym)/length(t);
		if m<0.5
			gc(i)=1;
		end

		if 0
			figure(1),clf,hold on
			B(31)=0;
			bar(B,C)
			plot(B,asym,'r')
			keyboard
		end

	end

end


A = saveAllResources(A);

