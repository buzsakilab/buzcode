function ConvertKilosort2Neurosuite(basepath,basename,rez)
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
% the folder under basepath with the data of interst.  Default = 1.
% basepath - directory path to the main recording folder with .dat and .xml
% as well as shank folders made by makeProbeMapKlusta2.m (default is
% current directory matlab is pointed to)
% basename - shared file name of .dat and .xml (default is last part of
% current directory path, ie most immediate folder name)
%
% Brendon Watson 2016

if ~exist('basepath','var')
    [~,basename] = fileparts(cd);
    basepath = cd;
end
if ~exist('rez','var');
    load(fullfile(basepath,'rez.mat'))
end


Nchan = rez.ops.Nchan;
connected    = ones(Nchan, 1);
xcoords      = ones(Nchan, 1);
ycoords      = (1:Nchan)';

par = LoadPar(fullfile(basepath,[basename '.xml']));
totalch = par.nChannels;
sbefore = 16;%samples before/after for spike extraction
safter = 24;%... could read from SpkGroups in xml
datpath = rez.ops.fbinary;

% [spikeTimes, ii] = sort(spikeTimes);
spktimes = uint64(rez.st3(:,1));
clu = uint32(rez.st3(:,2));
amplitudes = rez.st3(:,3);
pcFeatures = rez.cProjPC;
% pcFeatureInds = uint32(rez.iNeighPC);

%% do homework for assigning templates to shanks
% [~,shank]=fileparts(basepath);
templates = rez.Wraw;
% m = min(templates,[],2);%find the min value of each waveform on each channel
% [~,m] = min(m,[],1);%find which channel minimum is least
% m = squeeze(m);%which channel is minimum on each template
m = max(abs(templates),[],2);%find the most deviated value of each waveform on each channel
[~,m] = max(m,[],1);%find which channel has most deviated value
m = squeeze(m);%squeeze to 1d vector

grouplookup = rez.kcoords;
templateshankassignments = grouplookup(m);
allgroups = unique(grouplookup);
for groupidx = 1:length(allgroups)
    % for each group loop through, find all templates clus
    tgroup = allgroups(groupidx);%shank number
    ttemplates = find(templateshankassignments==tgroup);%which templates/clusters are in that shank
    tidx=ismember(clu,ttemplates);%find spikes indices in this shank
    tclu = clu(tidx);%extract template/cluster assignments of spikes on this shank
    tspktimes = spktimes(tidx);
    
    cluster_names = unique(tclu);
    t = find(rez.kcoords == tgroup);
    channellist = [];
    for ch = 1:length(par.SpkGrps)
        if ismember(t(1),par.SpkGrps(ch).Channels+1)
            channellist = par.SpkGrps(ch).Channels+1;
            continue
        end
    end
    if isempty(channellist)
        disp(['Cannot find spkgroup for group ' num2str(groupidx) ])
        continue
    end
    %% spike extraction from dat
    dat=memmapfile(datpath,'Format','int16');
    tsampsperwave = (sbefore+safter);
    ngroupchans = length(channellist);
    valsperwave = tsampsperwave * ngroupchans;
    wvforms_all=zeros(length(tspktimes)*tsampsperwave*ngroupchans,1,'int16');
    wvranges = zeros(length(tspktimes),ngroupchans);
    wvpowers = zeros(1,length(tspktimes));
    for j=1:length(tspktimes)
        try
            w = dat.data((double(tspktimes(j))-sbefore).*totalch+1:(double(tspktimes(j))+safter).*totalch);
            wvforms=reshape(w,totalch,[]);
            %select needed channels
            wvforms = wvforms(channellist,:);
    %         % detrend
    %         wvforms = floor(detrend(double(wvforms)));
            % median subtract
            wvforms = wvforms - repmat(median(wvforms')',1,sbefore+safter);
            wvforms = wvforms(:);
        catch
            disp(['Error extracting spike at sample ' int2str(double(tspktimes(j))) '. Saving as zeros']);
            disp(['Time range of that spike was: ' num2str(double(tspktimes(j))-sbefore) ' to ' num2str(double(tspktimes(j))+safter) ' samples'])
            wvforms = zeros(valsperwave,1);
        end

        %some processing for fet file
        wvaswv = reshape(wvforms,tsampsperwave,ngroupchans);
        wvranges(j,:) = range(wvaswv);
        wvpowers(j) = sum(sum(wvaswv.^2));

        lastpoint = tsampsperwave*ngroupchans*(j-1);
        wvforms_all(lastpoint+1 : lastpoint+valsperwave) = wvforms;
    %     wvforms_all(j,:,:)=int16(floor(detrend(double(wvforms)')));
        if rem(j,50000) == 0
            disp([num2str(j) ' out of ' num2str(length(tspktimes)) ' done'])
        end
    end
    clear dat
    wvranges = wvranges';

    %% Spike features
%     fets = h5read(tkwx,['/channel_groups/' num2str(shank) '/features_masks']);
%     fets = double(squeeze(fets(1,:,:)));
%     %mean activity per spike
%     fetmeans = mean(fets,1);
%     %find first pcs, make means of those... 
%     featuresperspike = 4;
%     firstpcslist = 1:featuresperspike:size(fets,1);
%     firstpcmeans = mean(fets(firstpcslist,:),1);
% 
%     nfets = size(fets,1)+1;
%     fets = cat(1,fets,fetmeans,firstpcmeans,wvpowers,wvranges,double(tspktimes'));
    fets = cat(1,wvpowers,wvranges,double(tspktimes'));
    fets = fets';
    % fets = cat(1,nfets,fets);

    %% writing to clu, res, fet, spk
    cluname = fullfile(basepath,[basename '.clu.' num2str(tgroup)]);
    resname = fullfile(basepath,[basename '.res.' num2str(tgroup)]);
    fetname = fullfile(basepath,[basename '.fet.' num2str(tgroup)]);
    spkname = fullfile(basepath,[basename '.spk.' num2str(tgroup)]);

    %clu
    % if ~exist(cluname,'file')
        tclu = cat(1,length(unique(tclu)),double(tclu));
        fid=fopen(cluname,'w'); 
    %     fprintf(fid,'%d\n',clu);
        fprintf(fid,'%.0f\n',tclu);
        fclose(fid);
        clear fid
    % end

    %res
    fid=fopen(resname,'w'); 
    % fprintf(fid,'%d\n',tspktimes);
    fprintf(fid,'%.0f\n',tspktimes);
    fclose(fid);
    clear fid

    %fet
    SaveFetIn(fetname,fets);

    %spk
    fid=fopen(spkname,'w'); 
    fwrite(fid,wvforms_all,'int16');
    fclose(fid);
    clear fid 

    disp(['Shank ' num2str(tgroup) ' done'])
end



function SaveFetIn(FileName, Fet, BufSize);

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

