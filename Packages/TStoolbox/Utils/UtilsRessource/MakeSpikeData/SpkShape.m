function A = truc(A)

A = getResource(A,'SpikeData');
A = getResource(A,'Tetrode');
A = registerResource(A,'SpkShape','tsdArray',{[],1},'spkShape','shape of spikes');


[p,ds,e] = fileparts(current_dir(A));

spkShape = {};

%Try to read spk & clu files for the 6 tetrodes

for i=1:6

	display(['Tetrode ' num2str(i)])

	try
		[sShape] = SpikeShape(ds,i,current_dir(A));
		if length(sShape)~= sum(TT==i)
			warning(['error between database and spk files'])
			display(lasterr)
			for j=1:sum(TT==i)
				spkShape = [spkShape;{tsd()}];
			end
		else
			spkShape = [spkShape;sShape];
		end
	catch
		warning(['No ' num2str(i) 'th TT on that day']);
		display(lasterr)
		for j=1:sum(TT==i)
			spkShape = [spkShape;{tsd()}];
		end
	end

end

spkShape = tsdArray(spkShape);

A = saveAllResources(A);

