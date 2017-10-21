function [fet clu spktimes wav]= ConvertKlusta2Matlab(shank,basepath,...
    basename,wvformExtract,saveFiles,numCores,keepMUA)
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
% keepMUA   - binary (0 or 1) to choose whether to keep MUA as a separate
%             cluster (#1)
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
% edited by Aza and Antonio 7/2017

% p = inputParser;
% addRequired(p,'shank',@isnumeric)
% addRequired(p,'basepath',@isstr)
% addRequired(p,'basename',@isstr)
%
% addParameter(p,'saveFiles',[],@isnumeric)
% addParameter(p,'numCores',[1],@isnumeric)
% p.parse(varargin{:})
%
% shank = p.Results.shank;
% basepath = p.Results.basepath;
% saveFiles = p.Results.saveFiles;
% numCores = p.Results.numCores;
% basename = p.Results.basename;


%
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

if saveFiles % make file names
    cluname = fullfile(basepath,[basename '.clu.' num2str(shank)]);
    resname = fullfile(basepath,[basename '.res.' num2str(shank)]);
    fetname = fullfile(basepath,[basename '.fet.' num2str(shank)]);
    spkname = fullfile(basepath,[basename '.spk.' num2str(shank)]);
end

% Start grabbing data
datpath = fullfile(basepath,[basename '.dat']);
tkwik = fullfile(basepath,num2str(shank),[basename '_sh' num2str(shank) '.kwik']);
tkwx = fullfile(basepath,num2str(shank),[basename '_sh' num2str(shank) '.kwx']);
kwikinfo = h5info(tkwik,['/channel_groups/' num2str(shank) '/clusters/main']);
clu = h5read(tkwik,['/channel_groups/' num2str(shank) '/spikes/clusters/main']);
cluster_names = unique(clu);

try % different klusta versions have different names for this field
    totalch = h5readatt(tkwik,'/application_data/spikedetekt','nchannels');
catch
    totalch = h5readatt(tkwik,'/application_data/spikedetekt','n_channels');
end
try % some files don't have these variables saved in them?
    sbefore = h5readatt(tkwik,'/application_data/spikedetekt','extract_s_before');
    safter = h5readatt(tkwik,'/application_data/spikedetekt','extract_s_after');
catch
    wvform_nsamples = h5readatt(tkwik,'/application_data/spikedetekt','waveforms_nsamples');
    sbefore = wvform_nsamples/2;
    safter = wvform_nsamples/2;
end
channellist = h5readatt(tkwik,['/channel_groups/' num2str(shank)],'channel_order')+1;

%% Getting spiketimes
spktimes = h5read(tkwik,['/channel_groups/' num2str(shank) '/spikes/time_samples']);

%% spike extraction from dat
wav = [];
tsampsperwave = (sbefore+safter);
ngroupchans = length(channellist);
valsperwave = tsampsperwave * ngroupchans;
if wvformExtract & ~exist(spkname)  % only create .spk file if it doesn't already exist
    disp('loading from .dat file')
    dat=memmapfile(datpath,'Format','int16');
    wvforms_all = zeros(length(spktimes),valsperwave,'int16');
    
    if numCores>1 % if using multiple cores, start parfor loop
        delete(gcp('nocreate')); parpool(numCores)
        parfor j=1:length(spktimes)
            try
                byteList = [];
                for i = 1:length(channellist)
                    ind1 = (double(spktimes(j)).*totalch); % correction for problems with high values
                    ind11 = int64(ind1) - int64(sbefore.*totalch) + 1;
                    ind2 = (double(spktimes(j)).*totalch);
                    ind22 = int64(ind2) + int64(safter.*totalch);
                    ind = ind11:ind22;
                    
                    byteList = [byteList; ind(double(channellist(i)):totalch:end)];
                end
                wvforms = dat.data(byteList);
            catch
                disp(['Error extracting spike at sample ' int2str(double(spktimes(j))) '. Saving as zeros']);
                disp(['Time range of that spike was: ' num2str(double(spktimes(j))-sbefore) ' to ' num2str(double(spktimes(j))+safter) ' samples'])
                wvforms = zeros(valsperwave,1,'int16');
            end
            if rem(j,50000) == 0
                disp([num2str(j) ' out of ' num2str(length(spktimes)) ' done'])
            end
            wvforms_all(j,:)=wvforms(:);
        end
    else % if we are using only one core, use a regular for loop
        for j=1:length(spktimes)
            try
                byteList = [];
                for i = 1:length(channellist)
                    ind1 = (double(spktimes(j)).*totalch); % correction for problems with high values
                    ind11 = int64(ind1) - int64(sbefore.*totalch) + 1;
                    ind2 = (double(spktimes(j)).*totalch);
                    ind22 = int64(ind2) + int64(safter.*totalch);
                    ind = ind11:ind22;
                    
                    byteList = [byteList; ind(double(channellist(i)):totalch:end)];
                end
                wvforms = dat.data(byteList);
            catch
                disp(['Error extracting spike at sample ' int2str(double(spktimes(j))) '. Saving as zeros']);
                disp(['Time range of that spike was: ' num2str(double(spktimes(j))-sbefore) ' to ' num2str(double(spktimes(j))+safter) ' samples'])
                wvforms = zeros(valsperwave,1,'int16');
            end
            if rem(j,50000) == 0
                disp([num2str(j) ' out of ' num2str(length(spktimes)) ' done'])
            end
            wvforms_all(j,:)=wvforms(:);
        end
        
    end
    
    
    clear dat
    wvforms_all = reshape(wvforms_all',size(wvforms_all,1)*size(wvforms_all,2),1);
    wav = reshape(wvforms_all,ngroupchans,tsampsperwave,[]);
    wav = permute(wav,[3 1 2]);
elseif wvformExtract & exist(spkname) % load from .spk if it exists
    disp('loading from .spk file')
    fid = fopen(spkname,'r');
    wav = fread(fid,[1 inf],'int16');
    wav = reshape(wav,ngroupchans,tsampsperwave,[]);
    wav = permute(wav,[3 1 2]);
end
%% Spike features

fets = h5read(tkwx,['/channel_groups/' num2str(shank) '/features_masks']);
fets = double(squeeze(fets(1,:,:)));
fet = fets';
fetMultiplier = double(intmax('int32'))/max(abs(fets(:))); % masked klustakwik has small floats that need to be expanded before rounding below
fet = fet .* fetMultiplier;

%% writing to clu, res, fet, spk
if saveFiles
    
    %clu
    for i=1:length(kwikinfo.Groups)  % sometimes .kwik has more groups than clu (empty?)..
        group(i) = h5readatt(tkwik,kwikinfo.Groups(i).Name,'cluster_group');
        temp = strsplit(kwikinfo.Groups(i).Name,'/');
        clusterID(i) = str2num(temp{length(temp)}); % fix for if channel ordering is mixed up (we don't know that the first clusters are noise/MUA)
    end
    
    cluster_names = unique(clu);
    for i=1:length(cluster_names)
        [a b] =ismember(cluster_names(i),double(clusterID));
        if isempty(b)
            warning('could not find cluster name?')
        end
        if keepMUA
            if group(b) == 0 || group(b) == 3
                clu(clu == clusterID(b)) = 0; % noise
            elseif group(b) == 1
                clu(clu == clusterID(b)) = 1; % MUA
            end
            
        else
            if group(b)~=2
                clu(clu == clusterID(b)) = 0;  % re-name unsorted and noise as MUA/Noise cluster for FMATToolbox
            end
        end
    end
    
    clu = cat(1,length(unique(clu)),double(clu));
    fid=fopen(cluname,'w');
    fprintf(fid,'%.0f\n',clu);
    fclose(fid);
    clear fid
    
    fid=fopen(resname,'w');
    fprintf(fid,'%.0f\n',spktimes);
    fclose(fid);
    clear fid
    
    SaveFetIn(fetname,fet);
    %spk
    
    if wvformExtract & ~exist(spkname)
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

