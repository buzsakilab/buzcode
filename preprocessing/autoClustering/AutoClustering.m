function AutoClustering(fbasename,elec,varargin)

% USAGE:
%     clu = AutoClustering(fbasename,elec)
% 
% AutoClustering automtically cleans the clu file defined by fbasename and
% electrode numbers. The program will look for the files
% fbasename.fet/res/spk.elec in the current folder.
% 
% INPUT:
%     fbasename: char array
%     elec: a vector of electrode numbers
%
% optional:
%     AutoClustering(fbasename,elec,dim)
%     where dim is the number of channels in electro group (if not
%     defined, will read the first line of the fet file
%     
% AutoClustering is meant to clean the output of KlustaKwik. The first
% thing it does is to separate electrical artifacts and MUA from putative
% isolated units. To do so, it sorts out units which have no clear
% refractory period (based on Hill, Mehta and Kleinfeld, J Neurosci.,
% 2012). Threshold can be set in the parameter section of this file
% ("Rogue spike threshold"). Then, it separates electrical
% artifats from MUA based on the assumption that electrical artifacts are
% highly correlated on thel different channels: the average waveform of at
% least one channel has to be different from the across-channel average
% waveform by a certrain amount of total variance (can be set in the
% parameter section, "Deviation from average spike threshold")
%
% Once the program has determined which of the clusters are putative
% isolated units, it tries to merge them based on waveform similarity
% (mahalanobis distance) and quality of the refractory period in the new
% merged cluster (or "Inter Common Spike Interval" from MS Fee et al. J 
% Neurosci. Meth., 1996)
%
% Adrien Peyrache, 2012
% David Tingley, 2017 updated for Intan system and new functionality

% The default behavior of this program is to take the output from KlustaKwik 
% or the newer klusta-3.0 (.kwik) algorithms, organize this data for easier 
% manual cluster cutting, and save to the klusta-3.0 (.kwik) format

if ~isempty(dir('*kwik'))
 if nargin < 1 
    basepath = pwd;
    name = dir('*alg*');
    s = split(name.name,'.');
    fbasename = s{1};
    elec = str2num(s{end});
 end
    
dbstop if error
% Parameters
% Number recording sites
if ~isempty(varargin)
    dim = varargin{1};
    dim = dim(:);
    if any(double(int16(dim))~=dim)
        error('Number of dimensions must be an integer')
    end
    if size(dim,1) ~= length(elec) && length(dim) ~=1
        error('Number of dimensions must be a vector of the same length as electrode vecotr or a single value')
    end
    if length(dim) == 1
        dim = dim*ones(length(elec),1);
    end
else
    dim = zeros(length(elec),1);
end

% an ugly list of heuristics that should eventually be cleaned up...

samplingRate = 20000;
% Load spike waveforms? Used for detection of electrical artifacts
loadspk = 0;
% Refractory period in msec
tR = 1.5./1000;
% Censored period in msec (specific to the spike detection process)
tC = 0.01./1000;
% Rogue spike threshold (for MUA); value between 0 an 1
rogThres = 1.5;
% Relative deviation (from 0 to 1) from average spike threshold (for electrical artifacts)
% =1000 => bypass it
devMinThres = .3;
devMaxThres = 10;

% other arbitrary thresholds..
powThresh = .85;
isoMinTresh = 1;
lratioMinThresh = 10e-5;
isoMaxTresh = 2000;
lratioMaxThresh = 1000;

% Do Merging?
doMerge = 1;
% overwrite clu file
rewriteclu= 1;
% Write a log file?
WriteLogFile = 1;

if WriteLogFile
    tic;
    log = [];
end

elec = elec(:)';
if length(elec)>1
    for eIx=1:length(elec)
        AutoClustering(fbasename,elec(eIx),dim(eIx));
    end
else
    % Load fet, clu, res, and waveforms
    fprintf('Sorting electrode %i of %s\n',elec,fbasename)
    if 1 %~exist([fbasename '_sh' num2str(elec) '.res.' num2str(elec)]) & ~exist([fbasename '.res.' num2str(elec)])
        if exist([fbasename '_sh' num2str(elec) '.kwik']) > 0
            tkwik = fullfile(pwd,[fbasename '_sh' num2str(elec) '.kwik']);
            cd ..; basepath = pwd; cd(num2str(elec));
        else
            tkwik = fullfile(pwd,num2str(elec),[fbasename '_sh' num2str(elec) '.kwik']);
            basepath = pwd; 
        end
        kwikinfo = h5info(tkwik,['/channel_groups/' num2str(elec) '/clusters/original']);
        if ~loadspk
            [fet clu res ~] = ConvertKlusta2Matlab(elec,basepath,fbasename,0,0,1);
        elseif loadspk 
           [fet clu res wav] = ConvertKlusta2Matlab(elec,basepath,fbasename,1,0,1);
        end
        clu_orig = clu;
        clu = double(clu);
        cluster_names = unique(clu);
    else
        error('check the below code for compatibility...')
        fet = LoadFeatures(fbasename,elec,dim);
        clu = load([fbasename '.clu.' num2str(elec)]);
        clu = clu(2:end);
        res = load([fbasename '.res.' num2str(elec)])/20;
        if loadspk
            wav = LoadSpikeWaveforms([fbasename '.spk.' num2str(elec)],dim,32);
        end
    end
    
    %Power of each waveform - UNUSED so far
%     pw = sqrt(sum(sum(wav.^2,2),3));
    % Percentile of the power
%     p = prctile(pw,0.99);
    nFeats = size(fet,2);  % to calculate accurate stats, we should have more rows than features for all cells
    
    %  devFromMean quantifies how much the spike are different on
    %  the different channels (what is the maximal distance from the averaged
    %  spike)

    for ii=1:length(cluster_names)
        rg = double(res(clu==cluster_names(ii)))./samplingRate;
        if loadspk
            m = squeeze(mean(wav(clu==cluster_names(ii),:,:)));
            y = pdist(m,'euclidean')./max(sqrt(sum(m.^2,2)));
            devFromMean(ii) = max(y);
        end
        fractRogue(ii) = FractionRogueSpk(rg,tR,tC);    
%         meanPw(ii) = nanmean(pw(clu==cluster_names(ii)))/nanmean(pw);   
    end
    if loadspk
        ff = find([devFromMean < devMinThres |...
             devFromMean > devMaxThres]);  % enforces variability across channels (muscle)
    else
        ff = [];
    end
    
    fff = find([fractRogue > rogThres]);  % enforces refractory period 

    for ii=1:length(cluster_names)
        if length(find(clu==cluster_names(ii))) > nFeats
            [L(ii) LRatio(ii)] = L_Ratio(fet,find(clu==cluster_names(ii)));
            iso(ii) = IsolationDistance(fet,find(clu==cluster_names(ii)));
        else
            L(ii)=nan;LRatio(ii)=nan;iso(ii)=nan; % not enough spikes to quantify anything...
        end
    end
    % Here let's toss out clusters that are clearly muscle/electrical
    % we have: fractRogue, devFromMean, LRatio, and iso to use
    f = find([LRatio<lratioMinThresh |...
            LRatio > lratioMaxThresh &...
            iso < isoMinTresh |...
            iso > isoMaxTresh]);          % finds large artefacts 
    noiseIx = find(ismember(clu,cluster_names(unique([f,ff,fff]))));

    
    % takes all muscle/electrical artefact and merges to a single
    % garbage cluster
    badCluIdx = unique(clu(noiseIx));
    for ii=1:length(badCluIdx)
        log = [log sprintf('%d -> %d; Looks like muscle or electrical artifact\n',badCluIdx(ii),0)];
    end
    clu(noiseIx) = 0;
    cluster_names = unique(clu);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Display intermediary results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if 0
        for i=1:length(cluster_names)
            subplot(2,2,1)
            m = squeeze(mean(wav(clu==cluster_names(i),:,:)));
            plot(m'); title(['Deviation from mean wf: ' num2str(devFromMean(i))])
            subplot(2,2,2)
            rg = double(res(clu==cluster_names(i)))./samplingRate;
            dt = diff(rg);
            hist(dt(dt<.1),1000); title(['Fraction Rogue: ' num2str(fractRogue(i))]);
            subplot(2,2,3)
            title(['L-Ratio: ' num2str(LRatio(i)) ' IsoDist: ' num2str(iso(i))])
            subplot(2,2,4)
            if ismember(cluster_names(i),unique(clu(noiseIx)))
                title('tossing');
            else
                title('keeping');
            end
%             title(num2str(d(:,i)'))
            pause
        end
    end
    
    %reorder clusters here
    [clu log] = renumberclu(clu,log);
 
    % Here we select only clusters that correspond to putative units and that
    % have at least 20 spikes (otherwise errormatrix calculation fails)
    h = hist(clu,length(unique(clu)))';
    goodCluIx = ismember(clu,find(cluster_names < 1000 & h>nFeats)); 
    goodCluIx(noiseIx) = 0;
     
    if length(unique(clu(goodCluIx)))>1  % if we have more than one putative cluster...
        newclu = clu(goodCluIx);
        newres = res(goodCluIx);
        newfet = fet(goodCluIx,:);
        if doMerge  % merge similar clusters which are neither noise nor MUA
            try
                [newclu mergehistory] = mergeclu_slow(newclu,newres,newfet,tR,tC,rogThres,samplingRate);
                % log changes
                log = [log sprintf('merge_slow.m was run\n')];
                for ii=1:size(mergehistory,1)
                    log = [log sprintf('%d + %d -> %d\n',mergehistory(ii,1),mergehistory(ii,2),mergehistory(ii,3))];
                end
            catch
                warning(['merging failed: ' lasterr])
            end            
        else   % if we're not going to try merging, the resort clusters by similarity
            em = errormatrix(fet(goodCluIx,:),newclu);
            ems = max(em,em');   
            y = squareform((1-ems)-diag(diag(1-ems)),'tovector');
            Z = linkage(y);
            [T,perm] = dendrogram_struct(Z,0);
            cluIx = unique(newclu);
            for ii=1:length(cluIx)
                newclu = updateclu(newclu,cluIx(perm(ii)),max(cluIx)+ii+1);
            end
        end
        clu(goodCluIx) = newclu; % re-insert any merges or resorting of cluster ID's        
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Display final results (requires function CrossCorr, not in the
    % toolbox)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if 0
    cluIx = unique(newclu(newclu>1));
    for x=1:length(cluIx)
        rgx = res(newclu==cluIx(x));
        [ax,b] =  CrossCorr(rgx*10,rgx*10,1,60);ax(b==0)=0;
        for y=x+1:length(cluIx)
            rgy = res(newclu==cluIx(y));
            [ay,b] =  CrossCorr(rgy*10,rgy*10,1,60);ay(b==0)=0;
            [fxy fcorr] = crossrefract(rgx,rgy);
            [h,b] = CrossCorr(rgx*10,rgy*10,1,60);ac(b==0)=0;
            figure(1),clf
            subplot(1,3,1)
                bar(b,ax,1)
             subplot(1,3,2)
                bar(b,h,1)
                title([fxy fcorr])
             subplot(1,3,3)
                bar(b,ay,1)
            pause
        end
    end
    end

     % reorder clusters here
    [clu log] = renumberclu(clu,log);
    disp(log)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Write new clu file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if rewriteclu
        kwikinfo = h5info(tkwik,['/channel_groups/' num2str(elec) '/clusters/main']);
        kwikinfo_original = h5info(tkwik,['/channel_groups/' num2str(elec) '/clusters/original']);
        
        if length(cluster_names) > length(kwikinfo_original.Groups)
%             kwikinfo = h5info(tkwik,['/channel_groups/' num2str(elec) '/clusters/main']);
            error('we added to the number of clusters?')
        end
        % rewrite cluster ID's
        h5write(tkwik,['/channel_groups/' num2str(elec) '/spikes/clusters/main'],uint32(clu));
        
%         I need to change autoclustering to selectively delete cluster groups that no longer have units....
        % delete all old cluster gropus from /main
%         fid = H5F.open(kwikinfo.Filename,'H5F_ACC_RDWR','H5P_DEFAULT');
%         for gg = 1:length(kwikinfo.Groups)
%             H5L.delete(fid,[kwikinfo.Groups(gg).Name],'H5P_DEFAULT')
%         end
%         H5F.close(fid)
        % copy cluster groups from /original to /main for new clusters
        fid = H5F.open(tkwik,'H5F_ACC_RDWR','H5P_DEFAULT');
        for i =1:length(cluster_names) % rewrite cluster groups
            try  % this may fail if something didn't get deleted... that's ok
            H5L.copy(fid,kwikinfo_original.Groups(i).Name,...
                            fid,[kwikinfo.Name '/' num2str(cluster_names(i))],...
                            'H5P_DEFAULT','H5P_DEFAULT')
            catch
                disp([ num2str(cluster_names(i)) ' cluster group already exists'])
            end
            if doMerge
                h5writeatt(tkwik,['/channel_groups/' num2str(elec) '/clusters/main/'...
                    num2str(cluster_names(i))],'cluster_group',3) % a human has not looked at these results, everything should be 'unsorted' (3 in .kwik format)
            end
        end
        H5F.close(fid)
        
        if exist([num2str(elec) '/nohup.out'])  && exist([num2str(elec) '/' fbasename '_sh' num2str(elec) '.klg.' num2str(elec)])
            if exist([fbasename '_sh' num2str(elec) '.kwik']) > 0
                fileID = fopen('nohup.out','a');
                klgID = fopen([fbasename '_sh' num2str(elec) '.klg.' num2str(elec)],'a');
            else
                fileID = fopen([num2str(elec) '/nohup.out'],'a');
                klgID = fopen([num2str(elec)  '/' fbasename '_sh' num2str(elec) '.klg.' num2str(elec)],'a');
            end
            fmt = 'this elec has been autoclustered';
            fprintf(fileID,fmt);
            fprintf(klgID,fmt);
            fclose(klgID);
            fclose(fileID); % write to both nohup and .klg. log files
        else
            klgID = fopen([fbasename '_sh' num2str(elec) '.klg.' num2str(elec)],'a');
            fmt = 'this elec has been autoclustered';
            fprintf(klgID,fmt);
            fclose(klgID);
            warning('could not find nohup.out log file')
        end
    end

    if WriteLogFile
        % Create (of overwrite) a log file
        fid = fopen([fbasename '.alg.' num2str(elec)],'w');
        time_tot = toc;
        fprintf(fid,[log 'That took %f seconds\n'],time_tot);
    end
end

end
