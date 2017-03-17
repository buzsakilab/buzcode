function [fet clu spktimes wav]= ConvertKlusta2Matlab(shank,basepath,basename,wvformExtract,saveFiles,numCores)
% USAGE
%
% [fet clu spktimes wav]= ConvertKlusta2Matlab(shank,basepath,basename,wvformExtract,saveFiles,numCores)
% 
% INPUTS
%
% shank - the shank number (as a number, not text) equalling the name of
%         the folder under basepath with the data of interst.  Default = 1.
% basepath - directory path to the main recording folder with .dat and .xml
%            as well as shank folders made by makeProbeMapKlusta2.m (default is
%            current directory matlab is pointed to)
% basename - shared file name of .dat and .xml (default is last part of
%            current directory path, ie most immediate folder name)
% wvformExtract - binary (0 or 1) to choose whether to extract raw waveforms
% saveFiles - binary (0 or 1) to choose whether to save clu, res, fet, wav(.spk) 
%             files or just return them as outputs
% numCores  - double between 1 and inf that sets the number of cores for a
%             parloop to run (set 2-12 on SSD or RAID0 systems, set to 1
%             for data stored on a single hard drive
%
%
% OUTPUTS
%
% fet - features of waveforms
% clu - cluster ID's for waveforms
% spktimes(res) - spike times of waveforms(in sample #'s, not seconds)
% wav- waveform shapes taken from raw .dat file
%
%Converts .kwik/kwx files from Klusta into klusters-compatible
% fet,res,clu,spk files.  Works on a single shank of a recording, assumes a
% 16bit .dat and an .xml file is present in "basepath" (home folder) and 
% that they are named basename.dat and basename.xml.  Also assumes that
% subdirectories in that basepath are made for each shank with names
% specified by numbers (ie 1,2,3,4..8).  In each shank folder should be
% .kwik and .kwx files made by klusta with names as follows:
% basename_sh[shankumber].kwik/kwx. 
%
%
% Brendon Watson 2016 
% edited by David Tingley 1/2017

% Handle inputs
if ~exist('shank','var')
    shank = 1;
end
if ~exist('basepath','var');
    basepath = cd;
end
if ~exist('saveFiles','var');
    saveFiles = 0;
end
if ~exist('numCores','var');
    numCores = 1;
end
if ~exist('basename','var');
    [~,basename] = fileparts(pwd);
elseif strcmp(basename, 'lfpfile') 
    d = dir('*lfp');
    basename = d.name(1:end-4);
end

% Start grabbing data
datpath = fullfile(basepath,[basename '.dat']);
tkwik = fullfile(basepath,num2str(shank),[basename '_sh' num2str(shank) '.kwik']);
tkwx = fullfile(basepath,num2str(shank),[basename '_sh' num2str(shank) '.kwx']);
kwikinfo = h5info(tkwik,['/channel_groups/' num2str(shank) '/clusters/main']);
clu = h5read(tkwik,['/channel_groups/' num2str(shank) '/spikes/clusters/main']);
cluster_names = unique(clu);

totalch = h5readatt(tkwik,'/application_data/spikedetekt','nchannels');
sbefore = h5readatt(tkwik,'/application_data/spikedetekt','extract_s_before');
safter = h5readatt(tkwik,'/application_data/spikedetekt','extract_s_after');
channellist = h5readatt(tkwik,['/channel_groups/' num2str(shank)],'channel_order')+1;

%% Getting spiketimes
spktimes = h5read(tkwik,['/channel_groups/' num2str(shank) '/spikes/time_samples']);

%% spike extraction from dat
wav = []; 
if wvformExtract
    dat=memmapfile(datpath,'Format','int16');
    tsampsperwave = (sbefore+safter);
    ngroupchans = length(channellist);
    valsperwave = tsampsperwave * ngroupchans;
    wvforms_all = zeros(length(spktimes),valsperwave,'int16');
    delete(gcp('nocreate')); parpool(numCores)
    parfor j=1:length(spktimes)
        try
            byteList = [];
            for i = 1:length(channellist)
                ind = (double(spktimes(j))-sbefore).*totalch+1:(double(spktimes(j))+safter).*totalch;
                byteList = [byteList; ind(double(channellist(i)):totalch:end)];
            end
            wvforms = dat.data(byteList);
        catch
            disp(['Error extracting spike at sample ' int2str(double(spktimes(j))) '. Saving as zeros']);
            disp(['Time range of that spike was: ' num2str(double(spktimes(j))-sbefore) ' to ' num2str(double(spktimes(j))+safter) ' samples'])
            wvforms = zeros(valsperwave,1);
        end
        if rem(j,50000) == 0
            disp([num2str(j) ' out of ' num2str(length(spktimes)) ' done'])
        end
        %some processing for fet file
        % wvaswv = reshape(wvforms,tsampsperwave,ngroupchans);
        % wvranges(j,:) = range(wvaswv);
        % wvpowers(j) = sum(sum(wvaswv.^2));
        % lastpoint = tsampsperwave*ngroupchans*(j-1);
        % wvforms_all(lastpoint+1 : lastpoint+valsperwave) = wvforms;
        wvforms_all(j,:)=wvforms(:);
    end
    clear dat
    wvforms_all = reshape(wvforms_all',size(wvforms_all,1)*size(wvforms_all,2),1);
    wav = reshape(wvforms_all,ngroupchans,tsampsperwave,[]);
    wav = permute(wav,[3 1 2]);
end
%% Spike features
fets = h5read(tkwx,['/channel_groups/' num2str(shank) '/features_masks']);
fets = double(squeeze(fets(1,:,:)));

%mean activity per spike
% fetmeans = mean(fets,1);
%find first pcs, make means of those... 
% featuresperspike = 4;
% firstpcslist = 1:featuresperspike:size(fets,1);
% firstpcmeans = mean(fets(firstpcslist,:),1);
% nfets = size(fets,1)+1;
% fets = cat(1,fets,fetmeans,firstpcmeans,wvpowers,wvranges,double(spktimes'));
% fets = cat(1,nfets,fets);

fet = fets';

%% writing to clu, res, fet, spk
if saveFiles
    cluname = fullfile(basepath,[basename '.clu.' num2str(shank)]);
    resname = fullfile(basepath,[basename '.res.' num2str(shank)]);
    fetname = fullfile(basepath,[basename '.fet.' num2str(shank)]);
    spkname = fullfile(basepath,[basename '.spk.' num2str(shank)]);

    %clu
    cluster_names = unique(clu);
    for i=1:length(cluster_names)
        group(i) = h5readatt(tkwik,kwikinfo.Groups(i).Name,'cluster_group');
        if group(i)~=2
        clu(clu == cluster_names(i)) = 0;  % re-name unsorted and noise as MUA/Noise cluster for FMATToolbox
        end
    end
    

    clu = cat(1,length(unique(clu)),double(clu));
    fid=fopen(cluname,'w'); 
    fprintf(fid,'%.0f\n',clu);
    fclose(fid);
    clear fid

    fid=fopen(resname,'w'); 
    % fprintf(fid,'%d\n',spktimes);
    fprintf(fid,'%.0f\n',spktimes);
    fclose(fid);
    clear fid
    
    SaveFetIn(fetname,fet);
    %spk
    if wvformExtract
        fid=fopen(spkname,'w'); 
        fwrite(fid,wvforms_all,'int16');
        fclose(fid);
        clear fid 
    end
    disp(['Shank ' num2str(shank) ' done'])
end

function SaveFetIn(FileName, Fet, BufSize)
if nargin<3 | isempty(BufSize)
    BufSize = inf;
end

nFeatures = size(Fet, 2);
formatstring = '%d';
for ii=2:nFeatures
  formatstring = [formatstring,'\t%d'];
end
formatstring = [formatstring,'\n'];

outputfile = fopen(FileName,'w');
fprintf(outputfile, '%d\n', nFeatures);

if isinf(BufSize)
    fprintf(outputfile,formatstring,round(Fet'));
else
    nBuf = floor(size(Fet,1)/BufSize)+1;
    
    for i=1:nBuf 
        BufInd = [(i-1)*nBuf+1:min(i*nBuf,size(Fet,1))];
        fprintf(outputfile,formatstring,round(Fet(BufInd,:)'));
    end
end

