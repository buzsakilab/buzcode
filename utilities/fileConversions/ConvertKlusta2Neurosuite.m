function ConvertKlusta2Neurosuite(shank,basepath,basename,numCores)
% Converts .kwik/kwx files from Klusta into klusters-compatible
% fet,res,clu,spk files.  Works on a single shank of a recording, assumes a
% 16bit .dat and an .xml file is present in "basepath" (home folder) and 
% that they are named basename.dat and basename.xml.  Also assumes that
% subdirectories in that basepath are made for each shank with names
% specified by numbers (ie 1,2,3,4..8).  In each shank folder should be
% .kwik and .kwx files made by klusta with names as follows:
% basename_sh[shankumber].kwik/kwx.  This 
% 
% Inputs:
% shank - the shank number (as a number, not text) equalling the name of
%           the folder under basepath with the data of interst.  Default = 1.
% basepath - directory path to the main recording folder with .dat and .xml
%           as well as shank folders made by makeProbeMapKlusta2.m (default is
%           current directory matlab is pointed to)
% basename - shared file name of .dat and .xml (default is last part of
%           current directory path, ie most immediate folder name)
% numCores  - double between 1 and inf that sets the number of cores for a
%             parloop to run (set 2-12 on SSD or RAID0 systems, set to 1
%             for data stored on a single hard drive

% Brendon Watson 2016
%
% Modified by Luke Sjulson, 2017-10

if ~exist('shank','var')
    shank = 1;
end
if ~exist('basepath','var');
    basepath = cd;
end
if ~exist('basename','var');
    [~,basename] = fileparts(cd);
elseif strcmp(basename, 'lfpfile')
    d = dir('*lfp');
    basename = d.name(1:end-4);
end
if ~exist('numCores','var');
    numCores = 1;
end


%% Start grabbing data
datpath = fullfile(basepath,[basename '.dat']);
tkwik = fullfile(basepath,num2str(shank),[basename '_sh' num2str(shank) '.kwik']);
tkwx = fullfile(basepath,num2str(shank),[basename '_sh' num2str(shank) '.kwx']);
% kwikinfo = h5info(tkwik,['/channel_groups/' num2str(shank) '/clusters/main']);
clu = h5read(tkwik,['/channel_groups/' num2str(shank) '/spikes/clusters/main']);
cluster_names = unique(clu);

try % different klusta versions have different names for this field
    totalch = h5readatt(tkwik,'/application_data/spikedetekt','nchannels');
catch
    totalch = h5readatt(tkwik,'/application_data/spikedetekt','n_channels');
end
sbefore = h5readatt(tkwik,'/application_data/spikedetekt','extract_s_before');
safter = h5readatt(tkwik,'/application_data/spikedetekt','extract_s_after');
channellist = h5readatt(tkwik,['/channel_groups/' num2str(shank)],'channel_order')+1;

%% find channels missing from spike groups becuase they were in "bad_channels.txt"
% this reflects what the xml and klusters expect but is not present
% will fill in data from missing channels with zeros
Par = LoadParameters(fullfile(basepath,[basename '.xml']));

spkgroupchannellist = Par.SpkGrps(shank).Channels + 1;
if length(spkgroupchannellist) > length(channellist)
    missingchannelidxs = find(spkgroupchannellist==setdiff(spkgroupchannellist,channellist));
else
    missingchannelidxs = [];
end

%% Getting spiketimes
spktimes = h5read(tkwik,['/channel_groups/' num2str(shank) '/spikes/time_samples']);
% spktimes = spktimes(1:10);

%% From Azahara: setting noise as 0 and MUA as 1
% when some files were corrupted I was able to still read them with this, not sure if all the corruption are like this one...
% Code that klustaviewa uses: 0 = noise, 1 = MUA, 2 = good
for ind = 1:length(cluster_names)
    kk=h5readatt(tkwik,['/channel_groups/' num2str(shank) '/clusters/main/' num2str(cluster_names(ind))],'cluster_group');
    cluster_group(ind) = kk(1,1);
    clear kk
    
    if cluster_group(ind) == 0 %NOISE
        clu(find(clu==cluster_names(ind))) = 0;
    elseif cluster_group(ind) == 1 %MUA
        clu(find(clu==cluster_names(ind))) = 1;
    end
end
clear ind cluster_group

%% spike extraction from dat
dat=memmapfile(datpath,'Format','int16');
tsampsperwave = (sbefore+safter);
ngroupchans = length(spkgroupchannellist);
valsperwave = tsampsperwave * ngroupchans;
wvforms_all = zeros(length(spktimes),valsperwave,'int16');
wvranges = zeros(length(spktimes),ngroupchans);
wvpowers = zeros(1,length(spktimes));

if numCores>1
    delete(gcp('nocreate')); parpool(numCores)
    parfor j=1:length(spktimes)
        % for j=1:length(spktimes)
        try
            byteList = [];
            for i = 1:length(channellist)
                ind = (double(spktimes(j))-sbefore).*totalch+1:(double(spktimes(j))+safter).*totalch;
                byteList = [byteList; ind(double(channellist(i)):totalch:end)];
            end
            wvforms = dat.data(byteList);
            %         w = dat.data((double(spktimes(j))-sbefore).*totalch+1:(double(spktimes(j))+safter).*totalch);
            %         wvforms=reshape(w,totalch,[]);
            %         %select needed channels
            %         wvforms = wvforms(channellist,:);
            % %         % detrend
            % %         wvforms = floor(detrend(double(wvforms)));
            for midx = 1:length(missingchannelidxs)%account for bad_channels that remained in spike groups in xml
                t = missingchannelidxs(midx);
                wvforms = cat(1,wvforms(1:t-1,:),zeros(1,size(wvforms,2)),wvforms(t:end,:));
            end
            % median subtract
            wvforms = wvforms - repmat(median(wvforms')',1,sbefore+safter);
        catch
            disp(['Error extracting spike at sample ' int2str(double(spktimes(j))) '. Saving as zeros']);
            disp(['Time range of that spike was: ' num2str(double(spktimes(j))-sbefore) ' to ' num2str(double(spktimes(j))+safter) ' samples'])
            wvforms = zeros(valsperwave,1);
        end
        wvforms_all(j,:)=wvforms(:);
        
        %some processing for fet file
        wvaswv = reshape(wvforms,ngroupchans,tsampsperwave)';
        wvranges(j,:) = range(wvaswv);
        wvpowers(j) = sum(sum(wvaswv.^2));
        
        %     lastpoint = tsampsperwave*ngroupchans*(j-1);
        wvforms_all(j,:)=wvforms(:);
        %     wvforms_all(lastpoint+1 : lastpoint+valsperwave) = wvforms;
        %     wvforms_all(j,:,:)=int16(floor(detrend(double(wvforms)')));
        if rem(j,50000) == 0
            disp([num2str(j) ' out of ' num2str(length(spktimes)) ' left'])
        end
    end
    
else % if using only one core, don't open a parpool
    for j=1:length(spktimes)
        % for j=1:length(spktimes)
        try
            byteList = [];
            for i = 1:length(channellist)
                ind = (double(spktimes(j))-sbefore).*totalch+1:(double(spktimes(j))+safter).*totalch;
                byteList = [byteList; ind(double(channellist(i)):totalch:end)];
            end
            wvforms = dat.data(byteList);
            %         w = dat.data((double(spktimes(j))-sbefore).*totalch+1:(double(spktimes(j))+safter).*totalch);
            %         wvforms=reshape(w,totalch,[]);
            %         %select needed channels
            %         wvforms = wvforms(channellist,:);
            % %         % detrend
            % %         wvforms = floor(detrend(double(wvforms)));
            for midx = 1:length(missingchannelidxs)%account for bad_channels that remained in spike groups in xml
                t = missingchannelidxs(midx);
                wvforms = cat(1,wvforms(1:t-1,:),zeros(1,size(wvforms,2)),wvforms(t:end,:));
            end
            % median subtract
            wvforms = wvforms - repmat(median(wvforms')',1,sbefore+safter);
        catch
            disp(['Error extracting spike at sample ' int2str(double(spktimes(j))) '. Saving as zeros']);
            disp(['Time range of that spike was: ' num2str(double(spktimes(j))-sbefore) ' to ' num2str(double(spktimes(j))+safter) ' samples'])
            wvforms = zeros(valsperwave,1);
        end
        wvforms_all(j,:)=wvforms(:);
        
        %some processing for fet file
        wvaswv = reshape(wvforms,ngroupchans,tsampsperwave)';
        wvranges(j,:) = range(wvaswv);
        wvpowers(j) = sum(sum(wvaswv.^2));
        
        %     lastpoint = tsampsperwave*ngroupchans*(j-1);
        wvforms_all(j,:)=wvforms(:);
        %     wvforms_all(lastpoint+1 : lastpoint+valsperwave) = wvforms;
        %     wvforms_all(j,:,:)=int16(floor(detrend(double(wvforms)')));
        if rem(j,50000) == 0
            disp([num2str(j) ' out of ' num2str(length(spktimes)) ' left'])
        end
    end
    
    
    
end


clear dat
wvforms_all = reshape(wvforms_all',size(wvforms_all,1)*size(wvforms_all,2),1);
wvranges = wvranges';

%% Spike features
fets = h5read(tkwx,['/channel_groups/' num2str(shank) '/features_masks']);
fets = double(squeeze(fets(1,:,:)));

% normalizing fets to be in a good numerical range for klusters (added by Luke, 2017-10)
f1 = size(fets,1);
f2 = size(fets,2);
fets = zscore(fets(:));
fets = reshape(fets, f1, f2);

fets(fets<-5) = -5; % trimming off everything >5 SDs from mean
fets(fets>5) = 5;

fets = fets.*1000;

% %mean activity per spike
% fetmeans = mean(fets,1);


%find first pcs, make means of those...
% featuresperspike = 4;%this is dumb, should read from file - not using it any more anyway.
% firstpcslist = 1:featuresperspike:size(fets,1);
% firstpcmeans = mean(fets(firstpcslist,:),1);

nfets = size(fets,1)+1;
% fets = cat(1,fets,fetmeans,firstpcmeans,wvpowers,wvranges,double(spktimes')); % these are a bunch of extra features that many people don't find helpful.
fets = cat(1,fets,double(spktimes')); % just the PCA projections and timestamps
fets = fets';
% fets = cat(1,nfets,fets);


%% writing to clu, res, fet, spk
cluname = fullfile(basepath,[basename '.clu.' num2str(shank)]);
resname = fullfile(basepath,[basename '.res.' num2str(shank)]);
fetname = fullfile(basepath,[basename '.fet.' num2str(shank)]);
spkname = fullfile(basepath,[basename '.spk.' num2str(shank)]);

%clu
clu = cat(1,length(unique(clu)),double(clu));
fid=fopen(cluname,'w');
fprintf(fid,'%.0f\n',clu);
fclose(fid);
clear fid

%res
fid=fopen(resname,'w');
fprintf(fid,'%.0f\n',spktimes);
fclose(fid);
clear fid

%fet
SaveFetIn(fetname,fets);

%spk
fid=fopen(spkname,'w');
fwrite(fid,wvforms_all,'int16');
fclose(fid);
clear fid

disp(['Shank ' num2str(shank) ' done'])


function SaveFetIn(FileName, Fet, BufSize);

if nargin<3 | isempty(BufSize)
    BufSize = inf;
end

nFeatures = size(Fet, 2);
formatstring = '%d';
for ii=2:nFeatures
    formatstring = [formatstring,'  %d'];
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

fclose(outputfile);