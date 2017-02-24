function S1 = SelectCells(S,minRate,epoch)

S1 = {};

for i=1:length(S)

	r = sum(Data(intervalCount(S{i},epoch)));
	if r>minRate
		S1 = [S1;{S{i}}];
	end

end

S1 = tsdArray(S1);