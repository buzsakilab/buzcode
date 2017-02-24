function ConvertSpykingCircus2Neurosuite(shank,basepath,basename)
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

if ~exist('shank','var')
    shank = 1;
end
if ~exist('basepath','var');
    [~,basename] = fileparts(cd);basepath = cd;
end

datpath = fullfile(basepath,[basename '.dat']);


% [~,shank]=fileparts(basepath);
tkwik = fullfile(basepath,num2str(shank),[basename '_sh' num2str(shank) '.kwik']);
tkwx = fullfile(basepath,num2str(shank),[basename '_sh' num2str(shank) '.kwx']);
clu = h5read(tkwik,['/channel_groups/' num2str(shank) '/spikes/clusters/main']);
cluster_names = unique(clu);

totalch = h5readatt(tkwik,'/application_data/spikedetekt','n_channels');
sbefore = h5readatt(tkwik,'/application_data/spikedetekt','extract_s_before');
safter = h5readatt(tkwik,'/application_data/spikedetekt','extract_s_after');
channellist = h5readatt(tkwik,['/channel_groups/' num2str(shank)],'channel_order')+1;

%% Getting spiketimes
spktimes = h5read(tkwik,['/channel_groups/' num2str(shank) '/spikes/time_samples']);
% spktimes = spktimes(1:10);

%% From Azahara: setting noise as 0 and MUA as 1
% Code that klustaviewa uses: 0 = noise, 1 = MUA, 2 = good
for ind = 1:length(cluster_names)
    %disp(cluster_names(ind));
    %when some files were corrupted I was able to still read them
    %with this, not sure if all the corruption are like this one...
    kk=h5readatt(tkwik,['/channel_groups/' num2str(shank) '/clusters/main/' num2str(cluster_names(ind))],'cluster_group');
    cluster_group(ind) = kk(1,1);
    clear kk

    if cluster_group(ind) == 0 %NOISE
        %spkts(find(clu==cluster_names(ind))) = []; %this option just if I want to remove the 0
        %waveform(:,:,find(clu==cluster_names(ind))) = [];
        %features(:,:,find(clu==cluster_names(ind))) = [];
        clu(find(clu==cluster_names(ind))) = 0;
    elseif cluster_group(ind) == 1 %MUA
        clu(find(clu==cluster_names(ind))) = 1;
    end
end
clear ind cluster_group 

%% spike extraction from dat
dat=memmapfile(datpath,'Format','int16');
tsampsperwave = (sbefore+safter);
ngroupchans = length(channellist);
valsperwave = tsampsperwave * ngroupchans;
wvforms_all=zeros(length(spktimes)*tsampsperwave*ngroupchans,1,'int16');
wvranges = zeros(length(spktimes),ngroupchans);
wvpowers = zeros(1,length(spktimes));
for j=1:length(spktimes)
    try
        w = dat.data((double(spktimes(j))-sbefore).*totalch+1:(double(spktimes(j))+safter).*totalch);
        wvforms=reshape(w,totalch,[]);
        %select needed channels
        wvforms = wvforms(channellist,:);
%         % detrend
%         wvforms = floor(detrend(double(wvforms)));
        % median subtract
        wvforms = wvforms - repmat(median(wvforms')',1,sbefore+safter);
        wvforms = wvforms(:);
    catch
        disp(['Error extracting spike at sample ' int2str(double(spktimes(j))) '. Saving as zeros']);
        disp(['Time range of that spike was: ' num2str(double(spktimes(j))-sbefore) ' to ' num2str(double(spktimes(j))+safter) ' samples'])
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
        disp([num2str(j) ' out of ' num2str(length(spktimes)) ' done'])
    end
end
clear dat
wvranges = wvranges';

%% Spike features
fets = h5read(tkwx,['/channel_groups/' num2str(shank) '/features_masks']);
fets = double(squeeze(fets(1,:,:)));
%mean activity per spike
fetmeans = mean(fets,1);
%find first pcs, make means of those... 
featuresperspike = 4;
firstpcslist = 1:featuresperspike:size(fets,1);
firstpcmeans = mean(fets(firstpcslist,:),1);

nfets = size(fets,1)+1;
fets = cat(1,fets,fetmeans,firstpcmeans,wvpowers,wvranges,double(spktimes'));
fets = fets';
% fets = cat(1,nfets,fets);


%% writing to clu, res, fet, spk
cluname = fullfile(basepath,[basename '.clu.' num2str(shank)]);
resname = fullfile(basepath,[basename '.res.' num2str(shank)]);
fetname = fullfile(basepath,[basename '.fet.' num2str(shank)]);
spkname = fullfile(basepath,[basename '.spk.' num2str(shank)]);

%clu
% if ~exist(cluname,'file')
    clu = cat(1,length(unique(clu)),double(clu));
    fid=fopen(cluname,'w'); 
%     fprintf(fid,'%d\n',clu);
    fprintf(fid,'%.0f\n',clu);
    fclose(fid);
    clear fid
% end

%res
fid=fopen(resname,'w'); 
% fprintf(fid,'%d\n',spktimes);
fprintf(fid,'%.0f\n',spktimes);
fclose(fid);
clear fid


%fet
SaveFetIn(fetname,fets);
% dlmwrite(fetname,nfets)
% dlmwrite(fetname,fets,'-append')
% % fid=fopen(fetname,'w'); 
% % fprintf(fid,'%i\n',nfets);
% % for a = 1:size(fets,1)
% %     fprintf(fid,'% i\n',fets(a,:));
% % end
% % fclose(fid);
% % clear fid


%spk
fid=fopen(spkname,'w'); 
fwrite(fid,wvforms_all,'int16');
fclose(fid);
clear fid 

disp(['Shank ' num2str(shank) ' done'])

% %--------------------------------------------------------------------------
% % Usage: ConvertPhy2KlustaSuite(filebase, shank_postfix)
% %
% % This script converts the result of the automated PHY clustering, making
% % it readable also by the older KlustaSuite's KlustaViewa manual clustering
% % application. It basically:
% % 1. appends the kwx and kwik file with some otherwise useless
% %    attribute entrys,
% % 2. copies some attributes to the old names which were renamed by the
% %    version change, and
% % 3. extracts the waveforms for waveform display from the dat file. Please
% %    note, that these waveforms are detrended raw waveforms and not
% %    filtered ones.
% % It is necessary to have the .kwik, .kwx, and the .dat file at the same
% % folder.
% % Shankpostfix is the postfix used for naming the kwik and kwx files for
% % individual shanks (e.g. '_sh5' in 'filebase_sh5.kwik'
% % WARNING: The script modifies the original kwx and kwik files, and makes
% % no backups. (Theoretically the PHY should still open them after the
% % modifications.)
% %
% % Antal Berenyi, 2015.11.20, drberenyi@gmail.com
% %--------------------------------------------------------------------------
% 
% 
% 
% %--------------------------------------------------------------
% % The script is using only low level H5F handling routines
% %--------------------------------------------------------------
% %Check if all file exists, quit if not.
% if ~exist('datpath','var')
%     datpath = ['../' filebase '.dat'];
% end
% 
% if (exist([filebase shank_postfix '.kwik'], 'file') ~= 2)
%     error([filebase shank_postfix '.kwik file is not accessible. Quitting.']);
% end
% if (exist([filebase shank_postfix '.kwx'], 'file') ~= 2) 
%     error([filebase shank_postfix '.kwx file is not accessible. Quitting.']);
% end
% if (exist(datpath, 'file') ~= 2)
%     error([filebase '.dat file is not accessible. Quitting.']);
% end
% 
% % Opening the kwik file
% disp(['Converting ' filebase shank_postfix ' from Phy format to KlustaSuite format...']);
% file=[filebase shank_postfix '.kwik'];
% fid=H5F.open(file,'H5F_ACC_RDWR','H5P_DEFAULT');
% 
% %--------------------------------------------------------------
% %--------------------------------------------------------------
% % The kwik file has to have some useless attributes at many places. Now we
% % collect these locations, and will write the attributes in a following
% % step in bulk.
% 
% places(1)={'/user_data'};
% places(2)={'/event_types'};
% %for all channel groups
% gid=H5G.open(fid,'/channel_groups');
% for i=0:H5G.get_num_objs(gid)-1;
%     places(end+1)={['/channel_groups/' H5G.get_objname_by_idx(gid,i) '/application_data']};
%     places(end+1)={['/channel_groups/' H5G.get_objname_by_idx(gid,i) '/user_data']};
%     gid2=H5G.open(fid,['/channel_groups/' H5G.get_objname_by_idx(gid,i) '/channels']);
%     for j=0:H5G.get_num_objs(gid2)-1;
%         places(end+1)={['/channel_groups/' H5G.get_objname_by_idx(gid,i) '/channels/' H5G.get_objname_by_idx(gid2,j) '/application_data']};
%         places(end+1)={['/channel_groups/' H5G.get_objname_by_idx(gid,i) '/channels/' H5G.get_objname_by_idx(gid2,j) '/application_data/klustaviewa']};
%         places(end+1)={['/channel_groups/' H5G.get_objname_by_idx(gid,i) '/channels/' H5G.get_objname_by_idx(gid2,j) '/application_data/spikedetekt']};
%         places(end+1)={['/channel_groups/' H5G.get_objname_by_idx(gid,i) '/channels/' H5G.get_objname_by_idx(gid2,j) '/user_data']};
%     end
%     H5G.close(gid2);
%     gid2=H5G.open(fid,['/channel_groups/' H5G.get_objname_by_idx(gid,i) '/cluster_groups/main']);
%     for j=0:H5G.get_num_objs(gid2)-1;
%         places(end+1)={['/channel_groups/' H5G.get_objname_by_idx(gid,i) '/cluster_groups/main/' H5G.get_objname_by_idx(gid2,j) '/user_data']};
%     end
%     H5G.close(gid2);
%     gid2=H5G.open(fid,['/channel_groups/' H5G.get_objname_by_idx(gid,i) '/cluster_groups/original']);
%     for j=0:H5G.get_num_objs(gid2)-1;
%         places(end+1)={['/channel_groups/' H5G.get_objname_by_idx(gid,i) '/cluster_groups/original/' H5G.get_objname_by_idx(gid2,j) '/user_data']};
%     end
%     H5G.close(gid2);
%     gid2=H5G.open(fid,['/channel_groups/' H5G.get_objname_by_idx(gid,i) '/clusters/main']);
%     for j=0:H5G.get_num_objs(gid2)-1;
%         places(end+1)={['/channel_groups/' H5G.get_objname_by_idx(gid,i) '/clusters/main/' H5G.get_objname_by_idx(gid2,j) '/user_data']};
%         places(end+1)={['/channel_groups/' H5G.get_objname_by_idx(gid,i) '/clusters/main/' H5G.get_objname_by_idx(gid2,j) '/quality_measures']};
%     end
%     H5G.close(gid2);
%     gid2=H5G.open(fid,['/channel_groups/' H5G.get_objname_by_idx(gid,i) '/clusters/original']);
%     for j=0:H5G.get_num_objs(gid2)-1;
%         places(end+1)={['/channel_groups/' H5G.get_objname_by_idx(gid,i) '/clusters/original/' H5G.get_objname_by_idx(gid2,j) '/user_data']};
%         places(end+1)={['/channel_groups/' H5G.get_objname_by_idx(gid,i) '/clusters/original/' H5G.get_objname_by_idx(gid2,j) '/quality_measures']};
%     end
%     H5G.close(gid2);
%     places(end+1)={['/channel_groups/' H5G.get_objname_by_idx(gid,i) '/spikes/features_masks']};
%     places(end+1)={['/channel_groups/' H5G.get_objname_by_idx(gid,i) '/spikes/waveforms_filtered']};
%     places(end+1)={['/channel_groups/' H5G.get_objname_by_idx(gid,i) '/spikes/waveforms_raw']};
%     objectname(i+1)={H5G.get_objname_by_idx(gid,i)};
% end
% gid=H5G.open(fid,'/recordings');
% for i=0:H5G.get_num_objs(gid)-1;
%     places(end+1)={['/recordings/' H5G.get_objname_by_idx(gid,i) '/user_data']};
%     places(end+1)={['/recordings/' H5G.get_objname_by_idx(gid,i) '/high']};
%     places(end+1)={['/recordings/' H5G.get_objname_by_idx(gid,i) '/low']};
%     objectnameRec(i+1)={H5G.get_objname_by_idx(gid,i)};
% end
% H5G.close(gid);
% H5F.close(fid);
% clear gid gid2 fid i j
% 
% %--------------------------------------------------------------
% %--------------------------------------------------------------
% % Write some completely useless, but needed attributes to thousands of
% % places in the kwik file. Places were collected before
% fid=H5F.open(file,'H5F_ACC_RDWR','H5P_DEFAULT');
% for j=1:length(places)
%     item=places{j};
%     entry={'CLASS'; 'TITLE'; 'VERSION'};
%     data={'GROUP'; ''; '1.0'};
%     gid=H5G.create(fid,item,0);
%     for i=1: length(data)
%         acpl_id = H5P.create('H5P_ATTRIBUTE_CREATE');
%         type_id = H5T.copy('H5T_C_S1');
%         H5T.set_size(type_id,max([1 length(data{i})]));
%         space_id = H5S.create('H5S_SCALAR');
%         attr_id = H5A.create(gid,entry{i},type_id,space_id,acpl_id);
%         if ~isempty(data{i})
%             H5A.write(attr_id,'H5ML_DEFAULT',data{i});
%         end
%         H5A.close(attr_id);
%     end
%     H5G.close(gid);
% end
% H5F.close(fid);
% 
% %--------------------------------------------------------------
% %--------------------------------------------------------------
% % Put links into the kwik file for the appropritate kwx file locations
% fid=H5F.open(file,'H5F_ACC_RDWR','H5P_DEFAULT');
% for j=1:length(objectname)
%     acpl_id = H5P.create('H5P_ATTRIBUTE_CREATE');
%     type_id = H5T.copy('H5T_C_S1');
%     data=['{kwx}/channel_groups/' objectname{j} '/waveforms_raw'];
%     H5T.set_size(type_id,max([1 length(data)]));
%     space_id = H5S.create('H5S_SCALAR');
%     gid=H5G.open(fid,['/channel_groups/' objectname{j} '/spikes/waveforms_raw']);
%     attr_id = H5A.create(gid,'hdf5_path',type_id,space_id,acpl_id);
%     H5A.write(attr_id,'H5ML_DEFAULT',data);
%     H5A.close(attr_id);
%     H5G.close(gid);
%     
%     data=['{kwx}/channel_groups/' objectname{j} '/waveforms_filtered'];
%     H5T.set_size(type_id,max([1 length(data)]));
%     space_id = H5S.create('H5S_SCALAR');
%     gid=H5G.open(fid,['/channel_groups/' objectname{j} '/spikes/waveforms_filtered']);
%     attr_id = H5A.create(gid,'hdf5_path',type_id,space_id,acpl_id);
%     H5A.write(attr_id,'H5ML_DEFAULT',data);
%     H5A.close(attr_id);
%     H5G.close(gid);
%     
%     data=['{kwx}/channel_groups/' objectname{j} '/features_masks'];
%     H5T.set_size(type_id,max([1 length(data)]));
%     space_id = H5S.create('H5S_SCALAR');
%     gid=H5G.open(fid,['/channel_groups/' objectname{j} '/spikes/features_masks']);
%     attr_id = H5A.create(gid,'hdf5_path',type_id,space_id,acpl_id);
%     H5A.write(attr_id,'H5ML_DEFAULT',data);
%     H5A.close(attr_id);
%     H5G.close(gid);
% end
% for j=1:length(objectnameRec) 
%     data=['{high.kwd}/recordings/' int2str(str2num(objectnameRec{j})) '/data'];
%     H5T.set_size(type_id,max([1 length(data)]));
%     space_id = H5S.create('H5S_SCALAR');
%     gid=H5G.open(fid,['/recordings/' int2str(str2num(objectnameRec{j})) '/high']);
%     attr_id = H5A.create(gid,'hdf5_path',type_id,space_id,acpl_id);
%     H5A.write(attr_id,'H5ML_DEFAULT',data);
%     H5A.close(attr_id);
%     H5G.close(gid);
%     
%     data=['{low.kwd}/recordings/' int2str(str2num(objectnameRec{j})) '/data'];
%     H5T.set_size(type_id,max([1 length(data)]));
%     space_id = H5S.create('H5S_SCALAR');
%     gid=H5G.open(fid,['/recordings/' int2str(str2num(objectnameRec{j})) '/low']);
%     attr_id = H5A.create(gid,'hdf5_path',type_id,space_id,acpl_id);
%     H5A.write(attr_id,'H5ML_DEFAULT',data);
%     H5A.close(attr_id);
%     H5G.close(gid);
%     
%     data=['{raw.kwd}/recordings/' int2str(str2num(objectnameRec{j})) '/data'];
%     H5T.set_size(type_id,max([1 length(data)]));
%     space_id = H5S.create('H5S_SCALAR');
%     gid=H5G.open(fid,['/recordings/' int2str(str2num(objectnameRec{j})) '/raw']);
%     attr_id = H5A.create(gid,'hdf5_path',type_id,space_id,acpl_id);
%     H5A.write(attr_id,'H5ML_DEFAULT',data);
%     H5A.close(attr_id);
%     H5G.close(gid);
% end
% H5F.close(fid);
% 
% %--------------------------------------------------------------
% %--------------------------------------------------------------
% % Write some paramteres of the experiment which are called differently in
% % phy
% fid=H5F.open(file,'H5F_ACC_RDWR','H5P_DEFAULT');
% acpl_id = H5P.create('H5P_ATTRIBUTE_CREATE');
% type_id = H5T.copy('H5T_NATIVE_LLONG');
% space_id = H5S.create('H5S_SCALAR');
% gid=H5G.open(fid,'/application_data/spikedetekt');
% attr_id = H5A.open(gid,'n_features_per_channel');
% data=H5A.read(attr_id);
% H5A.close(attr_id);
% attr_id = H5A.create(gid,'nfeatures_per_channel',type_id,space_id,acpl_id);
% H5A.write(attr_id,'H5ML_DEFAULT',data);
% H5A.close(attr_id);
% H5G.close(gid);
% 
% acpl_id = H5P.create('H5P_ATTRIBUTE_CREATE');
% type_id = H5T.copy('H5T_NATIVE_LLONG');
% space_id = H5S.create('H5S_SCALAR');
% gid=H5G.open(fid,'/application_data/spikedetekt');
% attr_id = H5A.open(gid,'extract_s_before');
% data=H5A.read(attr_id);
% H5A.close(attr_id);
% attr_id = H5A.open(gid,'extract_s_after');
% data=data+H5A.read(attr_id);
% H5A.close(attr_id);
% attr_id = H5A.create(gid,'waveforms_nsamples',type_id,space_id,acpl_id);
% H5A.write(attr_id,'H5ML_DEFAULT',data);
% H5A.close(attr_id);
% H5G.close(gid);
% 
% H5F.close(fid);
% clear fid gid i j data entry item space_id type_id acpl_id attr_id places
% 
% 
% %--------------------------------------------------------------
% %--------------------------------------------------------------
% % Create the matrices in the kwx file for the waveforms_raw and waveforms_filtered
% % do the kwx too
% file=[filebase shank_postfix '.kwx'];
% fid=H5F.open(file,'H5F_ACC_RDWR','H5P_DEFAULT');
% gid=H5G.open(fid,'/channel_groups');
% H5G.get_num_objs(gid);
% for i=1:H5G.get_num_objs(gid);
%     places(i)={['/channel_groups/' H5G.get_objname_by_idx(gid,i-1)]};
% end
% H5G.close(gid);
% 
% for i=1:length(places);
%     dset_id = H5D.open(fid,[places{i} '/features_masks']);
%     data = H5D.read(dset_id);
%     H5D.close(dset_id);
%     type_id = H5T.copy('H5T_NATIVE_SHORT');
%     h5_dims = [size(data,3) 32 size(data,2)/3];
%     h5_maxdims = h5_dims;
%     space_id = H5S.create_simple(3,h5_dims,h5_maxdims);
%     dcpl = 'H5P_DEFAULT';
%     dset_id = H5D.create(fid,[places{i} '/waveforms_raw'],type_id,space_id,dcpl);
%     H5D.close(dset_id);
%     dset_id = H5D.create(fid,[places{i} '/waveforms_filtered'],type_id,space_id,dcpl);
%     H5S.close(space_id);
%     H5T.close(type_id);
%     H5D.close(dset_id);
%     
% end
% H5F.close(fid);
% 
% %--------------------------------------------------------------
% %--------------------------------------------------------------
% %fill in the waveforms
% %read channel orders, spike times, and name of the shank from kwik file
% file=[filebase shank_postfix '.kwik'];
% fid=H5F.open(file,'H5F_ACC_RDWR','H5P_DEFAULT');
% for i=1:length(places);
%     gid=H5G.open(fid,places{i});
%     attr_id = H5A.open(gid,'channel_order');
%     chlist(i)={H5A.read(attr_id)};
%     H5A.close(attr_id);
%     H5G.close(gid);
%     dset_id = H5D.open(fid,[places{i} '/spikes/time_samples']);
%     spiketimes(i) = {H5D.read(dset_id)};
%     H5D.close(dset_id);
% end
% %get also total channel number, samples_before, samples_after
% gid=H5G.open(fid,'/application_data/spikedetekt');
% attr_id = H5A.open(gid,'n_channels');
% totalch=double(H5A.read(attr_id));
% H5A.close(attr_id);
% attr_id = H5A.open(gid,'extract_s_before');
% sbefore=double(H5A.read(attr_id));
% H5A.close(attr_id);
% attr_id = H5A.open(gid,'extract_s_after');
% safter=double(H5A.read(attr_id));
% H5A.close(attr_id);
% H5G.close(gid);
% H5F.close(fid);
% %use memmapfile .data iteratively to get all wavforms for all channels
% a=memmapfile(datpath,'Format','int16');
% for i=1:length(places);
%     wvforms_all(i)={zeros(length(spiketimes{i}), sbefore+safter ,length(chlist{i}),'int16')};
%     disp(['Extracting spikes waveforms for shank ' places{i} ' from dat file...']);
%     for j=1:length(spiketimes{i})
%         try
%         wvforms=reshape(a.data((double(spiketimes{i}(j))-sbefore).*totalch+1:(double(spiketimes{i}(j))+safter).*totalch),totalch,[]);
%         %select needed channels
%         wvforms=wvforms(chlist{i}+1,:);
%         %detrend them
%         wvforms_all{i}(j,:,:)=int16(floor(detrend(double(wvforms)')));
%         catch
%             disp(['Error by extracting spike at sample ' int2str(double(spiketimes{i}(j))) ', on shank ' places{i}]);
%         end
%     end
% end
% %put them in matrix in kwx file.
% file=[filebase shank_postfix '.kwx'];
% fid=H5F.open(file,'H5F_ACC_RDWR','H5P_DEFAULT');
% for i=1:length(places);
%     dset_id = H5D.open(fid,[places{i} '/waveforms_filtered']);
%     plist = 'H5P_DEFAULT';
%     H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',plist,permute(wvforms_all{i},[3 2 1]));
%     H5D.close(dset_id);
% end
% H5F.close(fid);
% clear chlist spiketimes totalch safter sbefore plist
% clear fid gid
% disp('Done.');
% 


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

