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
% highly correlated on the different channels: the average waveform of at
% least one channel has to be different from the across-channel average
% waveform by a certrain amount of total variance (can be set in the
% parameter section, "Deviation from average spike threshold")
%
%
% Once the program has determined which of the clusters are putative
% isolated units, it tries to merge them based on waveform similarity
% (mahalanobis distance) and quality of the refractory period in the new
% merged cluster (or "Inter Common Spike Interval" from MS Fee et al. J 
% Neurosci. Meth., 1996)
%
% Adrien Peyrache, 2012


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

% Load spike waveforms? Used for detection of electrical artifacts
loadspk = 1;
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

powThresh = .85;
% Use isolation dist and L ratio quality metrics to identify artifact
useMetrics = 0;

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
    if ~exist([fbasename '_sh' num2str(elec) '.res.' num2str(elec)]) & ~exist([fbasename '.res.' num2str(elec)])
        if exist([fbasename '_sh' num2str(elec) '.kwik']) > 0
            tkwik = fullfile(pwd,[fbasename '_sh' num2str(elec) '.kwik']);
            cd ..; basepath = pwd; cd(num2str(elec));
        else
            tkwik = fullfile(pwd,num2str(elec),[fbasename '_sh' num2str(elec) '.kwik']);
            basepath = pwd; 
        end
        kwikinfo = h5info(tkwik,['/channel_groups/' num2str(elec) '/clusters/main']);
        if ~loadspk
            [fet clu res ~] = ConvertKlusta2Neurosuite(elec,basepath,fbasename,0,0,1);
        elseif loadspk 
           [fet clu res wav] = ConvertKlusta2Neurosuite(elec,basepath,fbasename,1,0,1);
        end
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
    
    % number of spikes
    n = size(fet,1);

    %Power of each waveform - UNUSED so far
    pw = sqrt(sum(sum(wav.^2,2),3));
    % Percentile of the power
    p = prctile(pw,0.99);

    %  devFromMean quantifies how much the spike are different on
    %  the different channels (what is the maximal distance from the averaged
    %  spike)
    devFromMean = [];
    fractRogue = [];
    meanPw = [];

    for ii=1:length(cluster_names)
        rg = double(res(clu==cluster_names(ii)))./20000;
        if loadspk
            m = squeeze(mean(wav(clu==cluster_names(ii),:,:)));
            y = pdist(m,'euclidean')./max(sqrt(sum(m.^2,2)));
            devFromMean = [devFromMean,max(y)];
        end
        l = FractionRogueSpk(rg,tR,tC);    
        fractRogue = [fractRogue,l];
        meanPw = [meanPw,nanmean(pw(clu==cluster_names(ii)))/nanmean(pw)];   
    end
    f = [];
    ff = find([devFromMean < devMinThres |...
               devFromMean > devMaxThres]);  % enforces variability across channels (muscle)
    fff = find([fractRogue > rogThres]);  % enforces refractory period 

    if useMetrics
        for ii=1:length(cluster_names)
            if length(find(clu==cluster_names(ii))) > 50
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
    end
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
            rg = double(res(clu==cluster_names(i)))./20000;
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

    % Here we compute # of spike per cell. Some code for the errormatrix fails
    % when the cluster is defined by only a few samples. We'll put a
    % threshopld a bit later on the total # of spikes.
    h = hist(clu,double(unique(clu)));
    h = h(:);

    % Update the clu indices
%     if loadspk
%         for ii=1:length(noiseIx)
%             clu = updateclu(clu,noiseIx(ii),1000+ii);
%             log = [log sprintf('%d -> %d; Looks like electrical artifact\n',noiseIx(ii),1000+ii)];
%         end
%     end
%     for ii=1:length(muaIx)
%         clu = updateclu(clu,muaIx(ii),0+ii);
%         log = [log sprintf('%d -> %d; Looks like MUA\n',muaIx(ii),0+ii)];
%     end

    %reorder clusters here
    cluster_names = unique(clu);
    for i=1:length(cluster_names)
        if cluster_names(i) ~= 0
        clu(find(clu==cluster_names(i))) = i;
         log = [log sprintf('%d -> %d; reordering clusters\n',cluster_names(i),i)];
        end
    end
    
 
    % Here we select only clusters that correspond to putative units and that
    % have at least 20 spikes (otherwise errormatrix calculation fails)
    goodCluIx = ismember(clu,find(cluster_names < 1000 & h>20));
    goodCluIx(noiseIx) = 0;
    
    % if there is anything to compare...
    mergehistory = [];
 
    if length(unique(clu(goodCluIx)))>1
%         if max(clu(goodCluIx)) ~= max(clu(clu<1000))
%                 ix = clu == max(clu(goodCluIx));
%                 clu(ix) = max(clu(clu<1000))+1;
%         end
        newclu = clu(goodCluIx);
        newres = res(goodCluIx);
        newfet = fet(goodCluIx,:);
        if doMerge
            % merge similar clusters which are neither noise nor MUA
            try
                [newclu mergehistory] = mergeclu_slow(newclu,newres,newfet,tR,tC,rogThres);
                % log changes
                log = [log sprintf('merge_slow.m was run\n')];
                for ii=1:size(mergehistory,1)
                    log = [log sprintf('%d + %d -> %d\n',mergehistory(ii,1),mergehistory(ii,2),mergehistory(ii,3))];
                end
                if ~isempty(mergehistory) % for debugging
                   disp('clusters were merged') 
                end
            catch
                warning(['merging failed: ' lasterr])
            end            
        else
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
        clu(goodCluIx) = newclu;        
    end
    cluster_names = unique(clu);
    
    
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
    cluster_names = unique(clu);
    for i=1:length(cluster_names)
        if cluster_names(i) ~= 0
        clu(find(clu==cluster_names(i))) = i+2;
         log = [log sprintf('%d -> %d; reordering clusters\n',cluster_names(i),i)];
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Write new clu file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if rewriteclu
%         fid = fopen([fbasename '.clu.' num2str(elec)],'w');
%         fprintf(fid,'%i\n',[length(unique(clu));clu]);
%         fclose(fid);
        h5write(tkwik,['/channel_groups/' num2str(elec) '/spikes/clusters/main'],uint32(clu));
        if length(cluster_names) > length(kwikinfo.Groups)
            error('we added to the number of clusters?')
        end
        fid = H5F.open(tkwik,'H5F_ACC_RDWR','H5P_DEFAULT');
        for i =1:length(cluster_names) % rewrite cluster groups
            try
            H5L.copy(fid,[kwikinfo.Groups(i).Name],...
                fid,['/channel_groups/' num2str(elec) '/clusters/main/' num2str(cluster_names(i))],...
                'H5P_DEFAULT','H5P_DEFAULT')
            catch
                disp(['cluster ' num2str(i) ' already exists...'])
            end
        end
        H5F.close(fid)
        if exist([num2str(elec) '/nohup.out'])
            if exist([fbasename '_sh' num2str(elec) '.kwik']) > 0
                fileID = fopen('nohup.out','a');
            else
                fileID = fopen([num2str(elec) '/nohup.out'],'a');
            end
            fmt = 'this elec has been autoclustered';
            fprintf(fileID,fmt);
            fclose(fileID)
        else
            error('could not find nohup.out log file')
        end
    end

    if WriteLogFile
        % Create (of overwrite) a log file
        fid = fopen([fbasename '.alg.' num2str(elec)],'w');
        time_tot = toc;
        fprintf(fid,[log 'That took %f seconds\n'],time_tot);
    end
end
