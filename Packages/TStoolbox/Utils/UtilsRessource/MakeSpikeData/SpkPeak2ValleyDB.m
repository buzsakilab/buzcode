function A = truc(A)

A = getResource(A,'SpikeData');
A = getResource(A,'Tetrode');
A = registerResource(A,'P2v','numeric',{[],2},'p2v','peak to valley mean (1st column) and std (2nd) of spikes','mfile');

[p,ds,e] = fileparts(current_dir(A));

p2v = [];

%Try to read spk & clu files for the 6 tetrodes

for i=1:6

	display(['Tetrode ' num2str(i)])

	try
		t = SpikePeak2Valley(ds,i,current_dir(A));
		if size(t,1)~= sum(TT==i)
			warning(['error between database and spk files'])
			display(lasterr)
			p2v = [p2v;zeros(sum(TT==i),1)];
		else
			p2v = [p2v;t];
		end
	catch
		warning(['No ' num2str(i) 'th TT on that day']);
		display(lasterr)
		p2v = [p2v;zeros(sum(TT==i),1)];
	end

end

A = saveAllResources(A);

