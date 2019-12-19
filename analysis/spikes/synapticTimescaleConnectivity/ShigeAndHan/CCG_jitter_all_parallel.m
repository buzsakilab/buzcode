function ccgjitteroutput = CCG_jitter_all_parallel(basename,varargin)
% looks for connections between all cells in a given set of recordings with
% the input basename.  Calls Shige+Asohan's function ccg_jitter.m (modified
% by me.  Saves output to structure ccgjitteroutput.
% Brendon Watson 2012

%% gathering whether we'll plot, default is no plot, if user did not enter it.
if isempty(varargin);
    plotting = 0;
else
    plotting = varargin{1};
end

%% Read in basic cluster and spiking info, also recording info from hard drive
[spiket, spikeind, numclus, iEleClu, spikeph] = ReadEl4CCG2(basename);
Par = LoadPar([basename,'.xml']);

% 

%% Loop through each pair of cells, record output of shige's function into my struct
mtx = tril(meshgrid(1:numclus,1:numclus));
mtx(logical(eye(size(mtx))))=0;
[ixx,ixy]=find(mtx);


%% Preallocate
ccgjitteroutput(numclus,numclus).ccgR = []; 
ccgjitteroutput(numclus,numclus).tR = [];
ccgjitteroutput(numclus,numclus).GSPExc = [];
ccgjitteroutput(numclus,numclus).GSPInh = [];
ccgjitteroutput(numclus,numclus).ccgjMtx = [];
ccgjitteroutput(numclus,numclus).ccgjstats = [];
    

ccgjitteroutput(1,1).ConnectionsE = [];
ccgjitteroutput(1,1).ConnectionsI = [];
ccgjitteroutput(1,1).ECells = [];
ccgjitteroutput(1,1).ICells = [];
ccgjitteroutput(1,1).PossibleSameCells = [];
ccgjitteroutput(1,1).iEleClu = iEleClu;
ccgjitteroutput(1,1).SpikeIndices = spikeind;

%     iEleClu = iEleClu;
%     SpikeIndices = spikeind;

% warning off
% parfor ix = chosenindices;
%     a = ixx(ix);
%     b = ixy(ix);
%     [ccgR, tR, GSPExc, GSPInh, ccgjMtx, ccgjstats] = CCG_jitter(spiket,spikeind,a,b,Par.SampleRate,Par.SampleRate/1000,20,2,500,0.01,0);
% 
% %% Look for excitatory connections... will rule out incidents of same cell (delay=0), or other problems
%         findE = find(GSPExc); %find any time bins with significance for excitation from prior ccg_jitter.m
%         if ~isempty(findE) 
%             disp(['E Connection in pair ' num2str(a) ' and ' num2str(b)])
%             timeE = tR(findE);%recording times of positive bins
%             if ~isempty(find(timeE==0));
%                 PossibleSameCells(ix,:)=[a b];
%             end
%             
%             timeE = timeE(abs(timeE)<4 & abs(timeE)>1);%keep only those timestamps within synaptic window (1-4ms)
%             findE = findE(abs(timeE)<4 & abs(timeE)>1);%keep indices of above (clumsy, could change)
%             if ~isempty(findE); %if still any positives left (ie in time window)
%                 if timeE(1)>0;
%                     presyn=a; postsyn=b;
%                 elseif timeE(1)<0
%                     presyn=b; postsyn=a;
%                 end
%                 ConnectionsE(ix,:) = [presyn postsyn];%record as an exitatory cnxn
%                 ECells(ix) = presyn; %record this cell as excitatory
%             end
%         end
%         
% %% repeat with inhibitory cells
%         findI = find(GSPInh); %for inhib as well as excit pairings
%         if ~isempty(findI)
%             disp(['I Connection in pair ' num2str(a) ' and ' num2str(b)])
%             timeI = tR(findI);     
% %             if ~isempty(find(timeI==0));
% %                 ccgjitteroutput(1,1).PossibleSameCells(end+1,:)=[a b];
% %             end
%                 
%             timeI =timeI(abs(timeI)<4 & abs(timeI)>1);
%             findI = findI(abs(timeI)<4 & abs(timeI)>1);
%             if ~isempty(findI)
%                 if timeI(1)>0;
%                     presyn=a; postsyn=b;
%                 elseif timeI(1)<0
%                     presyn=b; postsyn=a;
%                 end
%                 ConnectionsI(ix,:) = [presyn postsyn];%record as an exitatory cnxn
%                 ICells(ix) = presyn;%record this cell as inhibitory
%             end
%         end
%         disp(['finished pair including ',num2str(a),' & ',num2str(b)])
% %     end
% end
% warning on

ixx = [1,3,5,7,1,2,5,6,1,2,3,4,1,2,3,4,1,2,3,1,2,1];
ixy = [2,4,6,7,3,4,7,8,4,5,6,7,5,6,7,8,6,7,8,7,8,8];

warning off

parfor ix = 1:length(ixx);
% chosenindices = 47:91;
% parfor ix = chosenindices;
    a = ixx(ix);
    b = ixy(ix);
    [ccgR{ix}, tR{ix}, GSPExc{ix}, GSPInh{ix}, ccgjMtx{ix}, ccgjstats{ix}] = CCG_jitter(spiket,spikeind,a,b,Par.SampleRate,Par.SampleRate/1000,20,2,500,0.01,0);
    disp(['Extracted CCGs for cells ',num2str(a),' & ',num2str(b)])
end
warning on

for ix = 1:length(ixx);
% for ix = chosenindices;
    a = ixx(ix);
    b = ixy(ix);
%     for b = (a+1):numclus
        ccgjitteroutput(a,b).ccgR =  ccgR{ix}; 
        ccgjitteroutput(a,b).tR =  tR{ix}; 
        ccgjitteroutput(a,b).GSPExc =  GSPExc{ix}; 
        ccgjitteroutput(a,b).GSPInh =  GSPInh{ix}; 
        ccgjitteroutput(a,b).ccgjMtx =  ccgjMtx{ix}; 
        ccgjitteroutput(a,b).ccgjstats =  ccgjstats{ix}; 
        
%% Look for excitatory connections... will rule out incidents of same cell (delay=0), or other problems
        findE = find(ccgjitteroutput(a,b).GSPExc); %find any time bins with significance for excitation from prior ccg_jitter.m
        if ~isempty(findE) 
            timeE = ccgjitteroutput(a,b).tR(findE);%recording times of positive bins
            if ~isempty(find(timeE==0));
                ccgjitteroutput(1,1).PossibleSameCells(end+1,:)=[a b];
            end
            
            timeE = timeE(abs(timeE)<4 & abs(timeE)>1);%keep only those timestamps within synaptic window (1-4ms)
            findE = findE(abs(timeE)<4 & abs(timeE)>1);%keep indices of above (clumsy, could change)
            if ~isempty(findE); %if still any positives left (ie in time window)
                if timeE(1)>0;
                    presyn=a; postsyn=b;
                elseif timeE(1)<0
                    presyn=b; postsyn=a;
                end
                ccgjitteroutput(1,1).ConnectionsE(end+1,:) = [presyn postsyn];%record as an exitatory cnxn
                ccgjitteroutput(1,1).ECells(end+1) = presyn; %record this cell as excitatory
            end
        end
        
%% repeat with inhibitory cells
        findI = find(ccgjitteroutput(a,b).GSPInh); %for inhib as well as excit pairings
        if ~isempty(findI)
            timeI = ccgjitteroutput(a,b).tR(findI);     
%             if ~isempty(find(timeI==0));
%                 ccgjitteroutput(1,1).PossibleSameCells(end+1,:)=[a b];
%             end
                
            timeI = timeI(abs(timeI)<4 & abs(timeI)>1);
            findI = findI(abs(timeI)<4 & abs(timeI)>1);
            if ~isempty(findI)
                if timeI(1)>0;
                    presyn=a; postsyn=b;
                elseif timeI(1)<0
                    presyn=b; postsyn=a;
                end
                ccgjitteroutput(1,1).ConnectionsI(end+1,:) = [presyn postsyn];%record as an exitatory cnxn
                ccgjitteroutput(1,1).ICells(end+1) = presyn;%record this cell as inhibitory
            end
        end
        disp(['finished pair including ',num2str(a),' & ',num2str(b)])
%     end
end



%% Gathering shank/cluster ID's over overlap cells (ie not just absolute cell number)
psc = ccgjitteroutput(1).PossibleSameCells;
if isempty(psc)
    ccgjitteroutput(1).PossibleSameCellIds = [];
else
    shanks = iEleClu(:,2);
    clus = iEleClu(:,3);
    ccgjitteroutput(1).PossibleSameCellIds = [shanks(psc(:,1)),clus(psc(:,1)),shanks(psc(:,2)),clus(psc(:,2))];
end

ccgjitteroutput(1,1).ECells = unique(ccgjitteroutput(1,1).ECells);
ccgjitteroutput(1,1).ICells = unique(ccgjitteroutput(1,1).ICells);


%% Classfying EE, EI, IE and II connections
ccgjitteroutput(1,1).ConnectionsEE = [];
ccgjitteroutput(1,1).ConnectionsEI = [];
ccgjitteroutput(1,1).ConnectionsEUnk = [];
ccgjitteroutput(1,1).ConnectionsIE = [];
ccgjitteroutput(1,1).ConnectionsII = [];
ccgjitteroutput(1,1).ConnectionsIUnk = [];

for a = 1:size(ccgjitteroutput(1,1).ConnectionsE)
    thisconnection = ccgjitteroutput(1,1).ConnectionsE(a,:);
    if ~isempty(find(ccgjitteroutput(1,1).ECells==thisconnection(2)));
        ccgjitteroutput(1,1).ConnectionsEE(end+1,:) = thisconnection;
    elseif ~isempty(find(ccgjitteroutput(1).ICells==thisconnection(2)));
        ccgjitteroutput(1,1).ConnectionsEI(end+1,:) = thisconnection;
    else
        ccgjitteroutput(1,1).ConnectionsEUnk(end+1,:) = thisconnection;
    end
end

for a = 1:size(ccgjitteroutput(1,1).ConnectionsI)
    thisconnection = ccgjitteroutput(1,1).ConnectionsI(a,:);
    if ~isempty(find(ccgjitteroutput(1,1).ECells==thisconnection(2)));
        ccgjitteroutput(1,1).ConnectionsIE(end+1,:) = thisconnection;
    elseif ~isempty(find(ccgjitteroutput(1).ICells==thisconnection(2)));
        ccgjitteroutput(1,1).ConnectionsII(end+1,:) = thisconnection;
    else
        ccgjitteroutput(1,1).ConnectionsIUnk(end+1,:) = thisconnection;
    end
end


%% plot outputs using other functions,  if the user chose do to so, default is no plot
if plotting
    try 
        CCG_jitter_plotpositives(ccgjitteroutput);
    catch
        if isempty(ccgjitteroutput(1).ConnectionsE) & isempty(ccgjitteroutput(1).ConnectionsI) 
            msgbox('There were no connections detected.')
        else
            msgbox('Unable to plot positive connections.  There were connections found, but there was an error.')
        end
    end

    try
        CCG_jitter_plot_possiblesamecells(ccgjitteroutput);
    catch
        if isempty(ccgjitteroutput(1).PossibleSameCells) 
            msgbox('There were no cells with significantly positive 0ms-lag CCGs indicating they may be the same.')
        else
            msgbox('Unable to plot possible same cells.  Such cells were found, but there was an error.')
        end
    end
end

%% save output to disk
%?smart)
save([basename '_ccgjitteroutput.mat'],'ccgjitteroutput');