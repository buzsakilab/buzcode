function ConvertMountainsort2Neurosuite(basepath,basename)
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

if ~exist('basepath','var');
    [~,basename] = fileparts(cd);
    basepath = cd;
end

datpath = fullfile(basepath,[basename '.dat']);
xmlpath = fullfile(basepath,[basename '.xml']);
originaldatapath = '/balrog_zpool/MountainSort/output/pipeAYA6Shank2--AYA6_day4/firings.mda';
datadir = fullfile(basepath);

% move some files to take care of issues with numeric names being rejected
% by h5 functions
% d = dir(fullfile(datadir,'*.result.hdf5'));
% oldresultspath = fullfile(datadir,d(end).name);
% newresultspath = fullfile(datadir,'result.result.hdf5');
% eval(['!cp ' oldresultspath ' ' newresultspath])

% d = dir(fullfile(datadir,'*.templates.hdf5'));
% oldtemplatespath = fullfile(datadir,d(end).name);
% newtemplatespath = fullfile(datadir,'templates.templates.hdf5');
% eval(['!cp ' oldtemplatespath ' ' newtemplatespath])

%% Get basic channel info
par = LoadXml(xmlpath);
totalch = par.nChannels;

sbefore = 16;%samples before/after for spike extraction
safter = 24;%... could read from SpkGroups in xml
if isfield(par.SpkGrps,'nSamples')
    if ~isempty(par.SpkGrps(1).nSamples);
        if isfield(par.SpkGrps,'PeakSample')
            if ~isempty(par.SpkGrps(1).PeakSample);
                sbefore = par.SpkGrps(1).PeakSample;
                safter = par.SpkGrps(1).nSamples - par.SpkGrps(1).PeakSample;
            end
        end
    end
end

grouplookup = zeros(totalch,1);
for a= 1:par.nElecGps
    grouplookup(par.ElecGp{a}+1) = a;
end
allgroups = unique(grouplookup);

%Grp 0 contain discared channels
allgroups(allgroups==0) = [];

%% get clusters and spike times
% clutimes = [];
% clugrps = [];
% for a = 1:numclus
%     eval(['clutimecell{' num2str(a) '}  = double(h5read(newresultspath,''/spiketimes/temp_' num2str(a-1) '''));'])
%     clutimes = cat(1,clutimes,clutimecell{a});
%     clugrps = cat(1,clugrps,a*ones(size(clutimecell{a})));
% end
% [spktimes,sidx] = sort(clutimes);
% clu = clugrps(sidx);

A = readmda(originaldatapath);
numclus = length(unique(A(3,:)));
spktimes = A(2,:);
clu = unique(A(3,:));

%% get templates
% templates_size = double(h5read(newtemplatespath, '/temp_shape'));
% N_e = templates_size(2);
% N_t = templates_size(1);
% temp_x = double(h5read(newtemplatespath, '/temp_x') + 1);
% temp_y = double(h5read(newtemplatespath, '/temp_y') + 1);
% temp_z = double(h5read(newtemplatespath, '/temp_data'));
% tmps = sparse(temp_x, temp_y, temp_z, templates_size(1)*templates_size(2), templates_size(3));
% templates_size = [templates_size(1) templates_size(2) templates_size(3)/2];
% for a = 1:numclus
%     templates(:,:,a) = full(reshape(tmps(:, a), templates_size(2), templates_size(1)));
% end

%% assign templates to shanks
% m = max(abs(templates),[],1);%find the most deviated value of each waveform on each channel
% [~,m] = max(m,[],2);%find which channel has most deviated value for each templnate
% m = squeeze(m);%squeeze to 1d vector
% 
% templateshankassignments = grouplookup(m);%for the list of maximal channels, which group is each in 

%% write shank-wise information
for groupidx = 1:length(allgroups)
    tgroup          = allgroups(groupidx);%shank number
    ttemplateidxs   = find(templateshankassignments==tgroup);%which templates/clusters are in that shank

    tidx            = ismember(clu,ttemplateidxs);%find spikes indices in this shank
    tclu            = clu(tidx);%extract template/cluster assignments of spikes on this shank
    tspktimes       = spktimes(tidx);
    
    %% 
    gidx            = find(grouplookup == tgroup);%find all channels in this group
    channellist     = [];
    for ch = 1:length(par.SpkGrps)
        if ismember(gidx(1),par.SpkGrps(ch).Channels+1)
            channellist = par.SpkGrps(ch).Channels+1;
            break
        end
    end
    if isempty(channellist)
        disp(['Cannot find spkgroup for group ' num2str(groupidx) ])
        continue
    end
    
    %% spike extraction from dat
    if groupidx == 1;
        dat = memmapfile(datpath,'Format','int16');
    end
    tsampsperwave   = (sbefore+safter);
    ngroupchans     = length(channellist);
    valsperwave     = tsampsperwave * ngroupchans;
    wvforms_all     = zeros(length(tspktimes)*tsampsperwave*ngroupchans,1,'int16');
%     wvranges        = zeros(length(tspktimes),ngroupchans);
%     wvpowers        = zeros(1,length(tspktimes));
    
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
%         wvranges(j,:) = range(wvaswv);
%         wvpowers(j) = sum(sum(wvaswv.^2));

        lastpoint = tsampsperwave*ngroupchans*(j-1);
        wvforms_all(lastpoint+1 : lastpoint+valsperwave) = wvforms;
    %     wvforms_all(j,:,:)=int16(floor(detrend(double(wvforms)')));
        if rem(j,100000) == 0
            disp([num2str(j) ' out of ' num2str(length(tspktimes)) ' done'])
        end
    end

    %% not doing fets for now
    
    %% writing to clu, res, fet, spk
    cluname = fullfile(basepath,[basename '.clu.' num2str(tgroup)]);
    resname = fullfile(basepath,[basename '.res.' num2str(tgroup)]);
%     fetname = fullfile(basepath,[basename '.fet.' num2str(tgroup)]);
    spkname = fullfile(basepath,[basename '.spk.' num2str(tgroup)]);
  %fet
%     SaveFetIn(fetname,fets);

    %clu
    tclu = cat(1,length(unique(tclu)),double(tclu));
    fid=fopen(cluname,'w'); 
    fprintf(fid,'%.0f\n',tclu);
    fclose(fid);
    clear fid

    %res
    fid=fopen(resname,'w'); 
    fprintf(fid,'%.0f\n',tspktimes);
    fclose(fid);
    clear fid

    %spk
    fid=fopen(spkname,'w'); 
    fwrite(fid,wvforms_all,'int16');
    fclose(fid);
    clear fid 

    disp(['Shank ' num2str(tgroup) ' done'])
    %end
    %end    
    
    
end
clear dat

%% Fets made the traditional way here
MakeClassicFet(basename,basepath)

%% Save original clus away, in case
mkdir(fullfile(basepath,'OriginalClus'))
copyfile(fullfile(basepath,[basename,'.clu.*']),fullfile(basepath,'OriginalClus'))

