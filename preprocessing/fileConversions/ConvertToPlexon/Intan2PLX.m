function  Intan2PLX(parameterSets,tet,varargin)


% function converts an Intan file to a plx file that can be read by Offline
% Sorter etc..http://www.plexon.com/software-downloads
%
% calls read_intan_data available at http://www.intantech.com/downloads.html
% requires spike detection software package from http://urut.ch/new/serendipity/index.php?/pages/osort.html
% calls scripts for writing PLX: checkPLXData, writePLXFile.
% readPLXHeaders is not called but is a useful reference
%
% INPUTS
%
% parameterSets:  %1 -> wavelet, 2->power
% tet: true -> groups electrodes by 4 , false ->treats electrode individually
%filename:
% ex: Intan2PLX(2,true,filename)





%define parameters of spike sort

params=[];
params.bandPass=[300 3000];
params.nrNoiseTraces=0;
params.prewhiten=0;
params.samplingFreq=40000; %intan sampling frequency
params.limit=100;

%set up filter for LFP
n = 4;
Wn =params.bandPass/(params.samplingFreq/2);
[b,a] = butter(n,Wn);
Hd=[];
Hd{1}=b;
Hd{2}=a;



switch( parameterSets )
    case 1
        %wavelet detection
        params.extractionThreshold = 0;
        params.alignMethod=2; %1 pos, 2 neg, 3 mixed
        params.peakAlignMethod=4; %1 findPeak, 2 none, 3 power, 4 mteo
        params.detectionMethod=5;
        dp3=[];
        dp3.scalesRange=[0.2 1.0];
        dp3.waveletName='bior1.5';
        params.detectionParams=dp3;
    case 2
        %power detection
        params.extractionThreshold = 6;
        params.alignMethod=2; %1 pos, 2 neg, 3 mixed
        params.peakAlignMethod=1; %1 findPeak, 2 none, 3 power, 4 mteo
        params.detectionMethod=1; %1 power, 5 WDM
        dp1=[];
        dp1.kernelSize=18;
        params.detectionParams=dp1;
    otherwise
        error('invalid parameter set nr');
end

if  length(varargin)>0
    fileList{1} = varargin{1};
else
    dirName=uigetdir;
    fileList = getAllIntanFiles(dirName);
    answer = inputdlg('Enter basefile');
    fileList = fileList(cellfun(@(a) any(regexp(a,answer{1})),fileList));
    
    %sort files by date
    for i=1:length(fileList)
    
    
    file = dir(fileList{i});
    filets(i,:)=datevec(file.date);
    end

    [a,b]=sortrows(filets);
    fileList=fileList(b);
end


%determine plexon file destination
[filename, pathname] = uiputfile('*.plx', 'Select a plx file');
fname = strcat(pathname, filename(1:end-4));

for ii = 1:length(fileList)
    %read intan data
    [t,amps,temp,aux] = read_intan_data(fileList{ii});
    data=zeros(length(temp),16); %hardcoded to 16 channels
    data(:,amps)=temp; %fill in dead channels as zeros
    data = data/1000;  %convert from uV to mV (plx standard)
    
    
    
    
    %%load spikes%%
    
    
    %initialize data output
    plexon.ts = nan(1e6,1);
    plexon.waves= nan(1e6,32); %approaching max memory limit
    plexon.units =nan(1e6,1);
    plexon.chans  =nan(1e6,1);
    
    
    
    %find spikes
    
    num_spks=0;
    
    %read 4 at a time if trodal data
    if tet
        skip = 4;
    else
        skip = 1;
    end
    
    
    for i=1:skip:size(data,2)
        
        
        %extract spikes from filtered LFP
        
        if tet
            spikeTimestamps_e{4}=[];
            filteredSignal_e{4}=[];
            for k=1:4
                [~, filteredSignal_e{k}, ~,~, spikeTimestamps_e{k}, ~, ~, ~] = extractSpikes(data(:,i+k-1), Hd, params );
                
                
            end
            
            %enslave signal from other channels to spikes on any electrode
            [spikeWaveforms,spikeTimestamps] = getTetrodeData(filteredSignal_e,spikeTimestamps_e,params.samplingFreq);
            
        else
            
            %not tested, likely will not work
            [~, filteredSignal, ~,spikeWaveforms, spikeTimestamps, ~, ~, ~] = extractSpikes(data(:,i), Hd, params );
        end
        
        %store vector array
        
        nspk=length(spikeTimestamps)*4;
        
        plexon.ts(num_spks+1:num_spks+nspk) = repmat(spikeTimestamps',4,1); %repeat ts for each unit
        plexon.waves(num_spks+1:num_spks+nspk,:) = spikeWaveforms; %all waveforms
        plexon.units(num_spks+1:num_spks+nspk) = zeros(nspk,1); %no unit assignments
        plexon.chans(num_spks+1:num_spks+nspk) = upSample(i:i+skip-1,length(spikeTimestamps)); %channel assignments
        
        
        num_spks=num_spks+nspk;
        
        
    end
    
    %remove excess variable space
    plexon.ts = plexon.ts(~isnan(plexon.ts));
    plexon.waves = plexon.waves(~any(isnan(plexon.waves),2),:);
    plexon.units = plexon.units(~isnan(plexon.units));
    plexon.chans = plexon.chans(~isnan(plexon.chans));
    
    
    %generate default headings
    headers = checkPLXData(plexon);
    
    %change headings specific to Intan
    headers.spikePreAmpGain=1; %signal already corrected for system gains
    headers.bitsperslowsample=16;
    headers.bitsperspikesample=16; %resolution of spike (y-axis)
    headers.ADFrequency=params.samplingFreq;
    
    headers.spikeMaxMagnitudeMV=max(abs(plexon.waves(:)));
    headers.slowMaxMagnitudeMV=max(abs(plexon.waves(:)));
    headers.waveformfreq=params.samplingFreq;
    plexon.waves = (plexon.waves * (.5*2^headers.bitsperslowsample))/ headers.spikeMaxMagnitudeMV; %convert to 16 bit
    plexon.waves(end-1,:)=plexon.waves(end,:);
    % write file
    
    writePLXFile([fname '_f' num2str(ii) '.plx'], plexon,headers)
    disp(['Saved in ' fname])
end

end




function [spikeWaveforms_tet,spikeTimestamps_tet] = getTetrodeData(filteredSignal,spikeTimestamps,fs)


filteredSignal=cell2mat(filteredSignal);
t_ind=1:size(filteredSignal,1);
spikes=false(size(t_ind));

%combine all spike times from a tetrode and sort
t=cell2mat(spikeTimestamps);
t=t(t>16 & t<size(filteredSignal,1)-17);
t = sort(t);

%find blocks with spikes on multiple channels or within 32 samples
if any(t)
epoch=repmat(-15:16,length(t),1)+repmat(t',1,32);
spikes(epoch(:))=true;
onSet=find(diff(spikes)==1);
offSet=find(diff(spikes)==-1);


%select biggest spike event in that block
for i=1:length(onSet)
    
    inT=t(t>onSet(i) & t<offSet(i));
    inT=unique(inT);
    inBlock = ismember(1:size(filteredSignal,1),inT);
    maxBlock = filteredSignal(inBlock,:);
    [r,~] = find(maxBlock==min(maxBlock(:)));
    ditch=inT(inT~=inT(r));
    
    t(ismember(t,ditch))=nan; %only keep biggest spike
    
end
t(isnan(t))=[];
t=unique(t);
spikeTimestamps_tet=t/fs; %convert to seconds
t=repmat(-15:16,length(t),1)+repmat(t',1,32);


temp1=filteredSignal(:,1); %e1
temp1=temp1(t);
temp2=filteredSignal(:,2); %e2
temp2=temp2(t);
temp3=filteredSignal(:,3); %e3
temp3=temp3(t);
temp4=filteredSignal(:,4); %e4
temp4=temp4(t);

spikeWaveforms_tet=[temp1;temp2;temp3;temp4];

else
    spikeWaveforms_tet=zeros(0,16);
    spikeTimestamps_tet=[];
end

end

function X= upSample (input, bins)
bins=round(bins);
input=input(:);
X=nan(length(input)*bins,1);


for i = 0:bins-1
    X(1+i:bins:end)=input;
end


end


