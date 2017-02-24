
function perf = performance(A,dset)

	for i=1:length(dset)
		A = getResource(A,'CorrectError',dset{i})
		ce = Data(correctError{1});
		perf(i) = sum(ce)/length(ce);
	end

end