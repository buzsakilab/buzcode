function [peak2peakWidth halfPeakWidth thirdPeakWidth] = SpikeWidth(day,TT,path)

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

peak2peakWidth = [];
halfPeakWidth = [];
thirdPeakWidth = [];
peakWidth = [];


for c=cellIdx

	cellSpk = waveforms(:,:,find(clu==c));
	
	for i=1:4	
		t = squeeze(cellSpk(i,:,:));
		v{i} = mean(t'); 
		[m(i),minPos(i)] = min(v{i}); %we want the biggest 2nd phase of spikes 
		baseL(i) = mean(v{i}([1,2,end-1,end]));
	end

	% on which tetrode is the 2nd peak the most negative?
	[dummy,gtt] = min(m);
	w = v{gtt};
	t = [1/fs:1/fs:nSamples/fs]*1000;

	wu = w - w(1);
	fs = fs;	

	% where is the positive peak?
	clear m maxPos
	[dummy,maxPos] = max(wu);

	%let's take the tail of spike shape on this tetrode
	waveTail = wu(maxPos:end);
	% where is the minimum?
	[m,minPos] = min(waveTail);

	if m<0
		%and the peak to peak widh is...
		peak2peakWidth = [peak2peakWidth;1000*minPos/fs];
		
		% Peak normalisation
		wu = waveTail/m;
		tu = t(maxPos:end);

		% Negative half peak width
		upHalfIx = find(wu>0.5);
		lowIx = upHalfIx(1);
		highIx = upHalfIx(end);

		if highIx<length(wu) & lowIx>1
	
			%   ...linear regression to find half crossing time
			t1 = tu(lowIx-1);
			t2 = tu(lowIx);
			y1 = wu(lowIx-1);
			y2 = wu(lowIx);
			tLow = (t2 - t1)/(y2 - y1)*(0.5 - y1) + t1; 
		
			t1 = tu(highIx);
			t2 = tu(highIx+1);
			y1 = wu(highIx);
			y2 = wu(highIx+1);
			tHigh = (t2 - t1)/(y2 - y1)*(0.5 - y1) + t1; 
	
			halfPeakWidth = [halfPeakWidth;tHigh-tLow];

		else
			halfPeakWidth = [halfPeakWidth;0];
		end

		% Negative 33% peak width
		upHalfIx = find(wu>0.34);
		lowIx = upHalfIx(1);
		highIx = upHalfIx(end);

		if highIx<length(wu) & lowIx>1

			%   ...linear regression to find half crossing time
			t1 = tu(lowIx-1);
			t2 = tu(lowIx);
			y1 = wu(lowIx-1);
			y2 = wu(lowIx);
			tLow = (t2 - t1)/(y2 - y1)*(0.5 - y1) + t1; 
		
			t1 = tu(highIx);
			t2 = tu(highIx+1);
			y1 = wu(highIx);
			y2 = wu(highIx+1);
			tHigh = (t2 - t1)/(y2 - y1)*(0.5 - y1) + t1; 
	
			thirdPeakWidth = [thirdPeakWidth;tHigh-tLow];

		else
			thirdPeakWidth = [thirdPeakWidth;0];
		end

	else
		peak2peakWidth = [peak2peakWidth;0];
		halfPeakWidth = [halfPeakWidth;0];
		thirdPeakWidth = [thirdPeakWidth;0];
	end

end