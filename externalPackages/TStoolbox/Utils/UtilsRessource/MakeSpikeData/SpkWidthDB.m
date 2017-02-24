function A = truc(A)

A = getResource(A,'SpikeData');
A = getResource(A,'Tetrode');
A = registerResource(A,'PpWidth','numeric',{[],1},'ppWidth','width of spikes','mfile');
A = registerResource(A,'HpWidth','numeric',{[],1},'hpWidth','width of spikes','mfile');
A = registerResource(A,'TpWidth','numeric',{[],1},'tpWidth','width of spikes','mfile');



[p,ds,e] = fileparts(current_dir(A));

ppWidth = [];
hpWidth = [];
tpWidth = [];

%Try to read spk & clu files for the 6 tetrodes

for i=1:6

	display(['Tetrode ' num2str(i)])

	try
		[t1 t2 t3] = SpikeWidth(ds,i,current_dir(A));
		if length(t1)~= sum(TT==i) | length(t2)~= sum(TT==i)
			warning(['error between database and spk files'])
			display(lasterr)
			ppWidth = [ppWidth;zeros(sum(TT==i),1)];
			hpWidth = [hpWidth;zeros(sum(TT==i),1)];
			tpWidth = [tpWidth;zeros(sum(TT==i),1)];
		else
			ppWidth = [ppWidth;t1];
			hpWidth = [hpWidth;t2];
			tpWidth = [tpWidth;t3];
		end
	catch
		warning(['No ' num2str(i) 'th TT on that day']);
		display(lasterr)
		ppWidth = [ppWidth;zeros(sum(TT==i),1)];
		hpWidth = [hpWidth;zeros(sum(TT==i),1)];
		tpWidth = [tpWidth;zeros(sum(TT==i),1)];
	end

end

A = saveAllResources(A);

