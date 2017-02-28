function ccgjitteroutput = CCG_jitter_classify(ccgjitteroutput,badcnxns)

if~exist('badcnxns','var')
    badcnxns = [];
end

if isfield(ccgjitteroutput,'ConnectionsE');
    ConnectionsE = ccgjitteroutput(1,1).ConnectionsE;
    ConnectionsI = ccgjitteroutput(1,1).ConnectionsI; 
    ECells = ccgjitteroutput(1,1).ECells;
    ICells = ccgjitteroutput(1,1).ICells;
    PossibleSameCells = ccgjitteroutput(1,1).PossibleSameCells;
%     SpikeIndices = ccgjitteroutput(1,1).SpikeIndices;
    ConnectionsEE = ccgjitteroutput(1,1).ConnectionsEE;
    ConnectionsEI = ccgjitteroutput(1,1).ConnectionsEI;
    ConnectionsEUnk = ccgjitteroutput(1,1).ConnectionsEUnk;
    ConnectionsIE = ccgjitteroutput(1,1).ConnectionsIE;
    ConnectionsII = ccgjitteroutput(1,1).ConnectionsII;
    ConnectionsIUnk = ccgjitteroutput(1,1).ConnectionsIUnk;
    temp = v2struct(ConnectionsE, ConnectionsI, ECells, ICells, PossibleSameCells,...
        ConnectionsEE, ConnectionsEI, ConnectionsEUnk,...
        ConnectionsIE,ConnectionsII,ConnectionsIUnk);
    ccgjitteroutput(1).oldConnectivity = temp;
%     ccgjitteroutput(1).oldConnectivity = v2struct(...
%         ccgjitteroutput(1,1).ConnectionsE, ccgjitteroutput(1,1).ConnectionsI,...
%         ccgjitteroutput(1,1).ECells, ccgjitteroutput(1,1).ICells,...
%         ccgjitteroutput(1,1).PossibleSameCells, ccgjitteroutput(1,1).SpikeIndices,...
%         ccgjitteroutput(1,1).ConnectionsEE, ccgjitteroutput(1,1).ConnectionsEI,...
%         ccgjitteroutput(1,1).ConnectionsEUnk, ccgjitteroutput(1,1).ConnectionsIE,...
%         ccgjitteroutput(1,1).ConnectionsII, ccgjitteroutput(1,1).ConnectionsIUnk);
    fields = {'ConnectionsE','ConnectionsI',...
        'ECells', 'ICells','PossibleSameCells',...
        'ConnectionsEE', 'ConnectionsEI','ConnectionsEUnk', 'ConnectionsIE',...
        'ConnectionsII', 'ConnectionsIUnk'};
    ccgjitteroutput=rmfield(ccgjitteroutput,fields);
    clear ConnectionsE ConnectionsI ECells ICells PossibleSameCells ...
        ConnectionsEE ConnectionsEI ConnectionsEUnk ConnectionsIE... 
        ConnectionsII ConnectionsIUnk 
end

ccgjitteroutput(1,1).ConnectionsE = [];
ccgjitteroutput(1,1).ConnectionsI = [];
ccgjitteroutput(1,1).ECells = [];
ccgjitteroutput(1,1).ICells = [];
ccgjitteroutput(1,1).PossibleSameCells = [];
% ccgjitteroutput(1,1).iEleClu = iEleClu;
% ccgjitteroutput(1,1).SpikeIndices = spikeind;

ccgjitteroutput(1,1).Connections0E = [];
ccgjitteroutput(1,1).Connections0I = [];

for a=1:size(ccgjitteroutput,2);      
    for b = (a+1):size(ccgjitteroutput,2);
        
%% Look for excitatory connections... will rule out incidents of same cell (delay=0), or other problems
        AllFoundEIndices = find(ccgjitteroutput(a,b).GSPExc); %find any time bins with significance for excitation from prior ccg_jitter.m
        if ~isempty(AllFoundEIndices) 
            AllFoundETimes = ccgjitteroutput(a,b).tR(AllFoundEIndices);%recording times of positive bins
            if ~isempty(find(AllFoundETimes==0));
                ccgjitteroutput(1,1).PossibleSameCells(end+1,:)=[a b];
            end
            
            GoodETimes = AllFoundETimes(abs(AllFoundETimes)<=4 & abs(AllFoundETimes)>=1);%keep only those timestamps within synaptic window (1-4ms)
%             GoodEIndices = AllFoundEIndices(abs(GoodETimes)<4 & abs(GoodETimes)>1);%keep indices of above (clumsy, could change)
            if ~isempty(GoodETimes); %if still any positives left (ie in time window)
                if GoodETimes(1)>0;
                    presyn=a; postsyn=b;
                elseif GoodETimes(1)<0
                    presyn=b; postsyn=a;
                end
                if isempty(find(ismember(badcnxns,[presyn postsyn],'rows')))
                    ccgjitteroutput(1,1).ConnectionsE(end+1,:) = [presyn postsyn];%record as an exitatory cnxn
                    ccgjitteroutput(1,1).ECells(end+1) = presyn; %record this cell as excitatory
                end
            end
            
            Good0ETimes = AllFoundETimes(abs(AllFoundETimes)==0);%keep only those timestamps at 0ms as a special category
%             find0E = find0E(abs(GoodETimes)==0);%keep indices of above (clumsy, could change)
            if ~isempty(Good0ETimes); %if still any positives left (ie in time window)
                if isempty(find(ismember(badcnxns,[a b],'rows'))) && isempty(find(ismember(badcnxns,[b a],'rows')))
                    ccgjitteroutput(1,1).Connections0E(end+1,:) = [a b];%record as an exitatory cnxn
                end
            end
        end
        
%% repeat with inhibitory cells
        AllFoundIIndices = find(ccgjitteroutput(a,b).GSPInh); %find any time bins with significance for excitation from prior ccg_jitter.m
        if ~isempty(AllFoundIIndices) 
            AllFoundITimes = ccgjitteroutput(a,b).tR(AllFoundIIndices);%recording times of positive bins
            
            GoodITimes = AllFoundITimes(abs(AllFoundITimes)<=4 & abs(AllFoundITimes)>=1);%keep only those timestamps within synaptic window (1-4ms)
%             GoodEIndices = AllFoundEIndices(abs(GoodETimes)<4 & abs(GoodETimes)>1);%keep indices of above (clumsy, could change)
            if ~isempty(GoodITimes); %if still any positives left (ie in time window)
                if GoodITimes(1)>0;
                    presyn=a; postsyn=b;
                elseif GoodITimes(1)<0
                    presyn=b; postsyn=a;
                end
                if isempty(find(ismember(badcnxns,[presyn postsyn],'rows')))
                    ccgjitteroutput(1,1).ConnectionsI(end+1,:) = [presyn postsyn];%record as an exitatory cnxn
                    ccgjitteroutput(1,1).ICells(end+1) = presyn; %record this cell as excitatory
                end
            end
%             
%             Good0ITimes = AllFoundITimes(abs(AllFoundITimes)==0);%keep only those timestamps at 0ms as a special category
% %             find0E = find0E(abs(GoodETimes)==0);%keep indices of above (clumsy, could change)
%             if ~isempty(Good0ITimes); %if still any positives left (ie in time window)
%                 ccgjitteroutput(1,1).Connections0I(end+1,:) = [a b];%record as an exitatory cnxn
%             endbasename

        end
        disp(b)
    end
end

%% Gathering shank/cluster ID's over overlap cells (ie not just absolute cell number)
psc = ccgjitteroutput(1).PossibleSameCells;
if isempty(psc)
    ccgjitteroutput(1).PossibleSameCellIds = [];
else
%    numclus shanks = iEleClu(:,2);
%     clus = iEleClu(:,3);
%     ccgjitteroutput(1).PossibleSameCellIds = [shanks(psc(:,1)),clus(psc(:,1)),shanks(psc(:,2)),clus(psc(:,2))];
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