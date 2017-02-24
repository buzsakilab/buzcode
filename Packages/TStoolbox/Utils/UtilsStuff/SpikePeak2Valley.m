function p2v = SpikeWidth(day,TT,path)

%Parameters

nChannels = 4;
nSamples = 34;
fs = 25000;


spkFile = [path filesep day '.spk.' num2str(TT)];
parFile = [path filesep day '.par.' num2str(TT)];
cluFile = [path filesep day '.clu.' num2str(TT)];

% Here, open parFile to get number of recorded points per spikes.

try
	par = fopen(parFile);
	for i=1:5
		l = fgetl(par);
	end
	nSamples = str2num(l(1:2));

catch
	warning('file par not present')
end

%try to load clu file before unzipping spk file. If doen't not exist generates an error

if exist(cluFile,'file')
	clu = load(cluFile);
else
	error([cluFile ' doesn''t exist']);
end

if exist([spkFile '.gz'])
	display(['unzipping file ' spkFile]);
	eval(['!gunzip ' spkFile '.gz']);
elseif exist(spkFile)
	display(['already uncompressed ' spkFile]);
else
	error([spkFile ' doesn''t exist']);
end

%first element of clu is # of clusters
cellIdx = 2:(clu(1)-1);
clu = clu(2:end);


% Load spike waveforms
try
	waveforms = LoadSpikeWaveforms(spkFile,nChannels,nSamples);
catch
	nSamples = 32;
	waveforms = LoadSpikeWaveforms(spkFile,nChannels,nSamples);
end
waveforms(:,:,find(clu<2)) = [];
clu(find(clu<2))=[];

display(['zipping file ' spkFile]);
eval(['!gzip ' spkFile]);

p2v = [];

for c=cellIdx

	cellSpk = waveforms(:,:,find(clu==c));
	
	for i=1:4	
		t = double(squeeze(cellSpk(i,:,:)));
		m = mean(t'); 

		[maxVal(i),maxPos(i)] = max(m);
		waveTail = m(maxPos(i):end);
		[minVal(i),minPos(i)] = min(waveTail);
	end

	% on which tetrode is the P2V the biggest?
	
	[p2vMean p2vIx] = max(maxVal-minVal);

	t = double(squeeze(cellSpk(p2vIx,:,:)))';
	mx = maxPos(p2vIx);
	mn = minPos(p2vIx)+maxPos(p2vIx)-1;
	p = (t(:,mx)-t(:,mn));
	p2vStd = std(p);

	p2v = [p2v;[p2vMean p2vStd]];

end