function funcsynapses = FindSynapseToStruct(funcsynapses,badcnxns,widecnxns)
% Identifies and stores in a structure types of connections and cell types
% based on connectivity, based on Make_FindSynapse_bw. Takes out bad
% connections and also re-categorizes wide connections.
%
% "funcsynapses" must already have in it the following:
% funcsynapses.CnxnTimes = synTimes;
% funcsynapses.CnxnWeightsZ = synWeightsZ;
% funcsynapses.CnxnWeightsR = synWeightsR;
% funcsynapses.PairUpperThreshs = triu(synSig);
% funcsynapses.PairLowerThreshs = tril(synSig);

if~exist('badcnxns','var') && ~exist('widecnxns','var')
    badcnxns = [];
    widecnxns = [];
else %ie if badcnxns was an input...gonna store away old connectivity metadata
    ConnectionsE = funcsynapses.ConnectionsE;
    ConnectionsI = funcsynapses.ConnectionsI; 
    ECells = funcsynapses.ECells;
    ICells = funcsynapses.ICells;
    ConnectionsEE = funcsynapses.ConnectionsEE;
    ConnectionsEI = funcsynapses.ConnectionsEI;
    ConnectionsEUnk = funcsynapses.ConnectionsEUnk;
    ConnectionsIE = funcsynapses.ConnectionsIE;
    ConnectionsII = funcsynapses.ConnectionsII;
    ConnectionsIUnk = funcsynapses.ConnectionsIUnk;

    if isfield(funcsynapses,'BadConnections');
        BadConnections = funcsynapses.BadConnections;
    else
        BadConnections = [];
    end
    if isfield(funcsynapses,'WideConnections');
        WideConnections = funcsynapses.WideConnections;
    else
        WideConnections = [];
    end
    
    temp = v2struct(ConnectionsE, ConnectionsI, ECells, ICells,...
        ConnectionsEE, ConnectionsEI, ConnectionsEUnk,...
        ConnectionsIE,ConnectionsII,ConnectionsIUnk,BadConnections,WideConnections);
    if ~isfield(funcsynapses,'OriginalConnectivity')
        newfieldname = 'OriginalConnectivity';
    else
        newfieldname = 'MostRecentConnectivity';
    end
    eval(['funcsynapses.' newfieldname '=temp;'])
%     fields = {'ConnectionsE','ConnectionsI','ECells', 'ICells',...
%         'ConnectionsEE', 'ConnectionsEI','ConnectionsEUnk', 'ConnectionsIE',...
%         'ConnectionsII', 'ConnectionsIUnk', 'BadConnections'};
%     funcsynapses=rmfield(funcsynapses,fields);
    clear ConnectionsE ConnectionsI ECells ICells...
        ConnectionsEE ConnectionsEI ConnectionsEUnk ConnectionsIE... 
        ConnectionsII ConnectionsIUnk BadConnections
end


funcsynapses.ConnectionsE = [];
funcsynapses.ConnectionsI = [];
funcsynapses.ECells = [];
funcsynapses.ICells = [];
funcsynapses.Connections0E = [];
funcsynapses.Connections0I = [];
funcsynapses.ConnectionsEE = [];
funcsynapses.ConnectionsEI = [];
funcsynapses.ConnectionsEUnk = [];
funcsynapses.ConnectionsIE = [];
funcsynapses.ConnectionsII = [];
funcsynapses.ConnectionsIUnk = [];
funcsynapses.ConnectionsEI = [];
funcsynapses.ConnectionsEUnk = [];
funcsynapses.ConnectionsIE = [];
funcsynapses.ConnectionsII = [];
funcsynapses.ConnectionsIUnk = [];
funcsynapses.BadConnections = badcnxns;
funcsynapses.WideConnections = widecnxns;

%% Find E connections (commented out was for finding correct pre/post direction, but Adrien already did it
% [refidx,nrefidx] = find(synWeightsZ>0);%grabbing excitatory cnxns
% reversed = synTimes(synWeightsZ>0)<0;
% for a = 1:size(refidx,1);
%     if reversed(a)
%         pre = nrefidx(a);
%         post = refidx(a);
%     else
%         pre = refidx(a);
%         post = nrefidx(a);
%     end
%     funcsynapses.ConnectionsE(a,:) = [pre post];
% end
[pre,post] = find(funcsynapses.CnxnWeightsZ>0);%grabbing excitatory cnxns
if ~isempty(badcnxns)
    badidx = [];
    for a = 1:size(pre,1);
        if ismember([pre(a) post(a)],badcnxns,'rows')
            badidx(end+1) = a;
        end
    end
    pre(badidx) = [];
    post(badidx) = [];
end
% No longer eliminate wide... will be indep from bad now
% if ~isempty(widecnxns)
%     wideidx = [];
%     for a = 1:size(pre,1);
%         if ismember([pre(a) post(a)],widecnxns,'rows')
%             wideidx(end+1) = a;
%         end
%     end
%     pre(wideidx) = [];
%     post(wideidx) = [];
% end

funcsynapses.ConnectionsE = [pre post];

%% Find I connections (commented out was for finding correct pre/post direction, but Adrien already did it
% [refidx,nrefidx] = find(synWeightsZ<0);%grabbing inhibitory cnxns
% reversed = synTimes(synWeightsZ<0)<0;
% for a = 1:size(refidx,1);
%     if reversed(a)
%         pre = nrefidx;
%         post = refidx;
%     else
%         pre = refidx;
%         post = nrefidx;
%     end
%     funcsynapses.ConnectionsI = [pre post];
% end
[pre,post] = find(funcsynapses.CnxnWeightsZ<0);%grabbing inhibitory cnxns
if ~isempty(badcnxns)
    badidx = [];
    for a = 1:size(pre,1);
        if ismember([pre(a) post(a)],badcnxns,'rows')
            badidx(end+1) = a;
        end
    end
    pre(badidx) = [];
    post(badidx) = [];
end
% if ~isempty(widecnxns)
%     wideidx = [];
%     for a = 1:size(pre,1);
%         if ismember([pre(a) post(a)],widecnxns,'rows')
%             wideidx(end+1) = a;
%         end
%     end
%     pre(wideidx) = [];
%     post(wideidx) = [];
% end
funcsynapses.ConnectionsI = [pre post];


%% Simple ID'ing of each cell
if ~isempty(funcsynapses.ConnectionsE)
    funcsynapses.ECells = unique(funcsynapses.ConnectionsE(:,1));
end

if ~isempty(funcsynapses.ConnectionsI)
    funcsynapses.ICells = unique(funcsynapses.ConnectionsI(:,1));
end

funcsynapses.EIProblemCells = intersect(funcsynapses.ECells,funcsynapses.ICells);

%% Classfying EE, EI, IE and II connections
if ~isempty(funcsynapses.ConnectionsE)
    for a = 1:size(funcsynapses.ConnectionsE)
        thisconnection = funcsynapses.ConnectionsE(a,:);
        if ~isempty(find(funcsynapses.ECells==thisconnection(2)));
            funcsynapses.ConnectionsEE(end+1,:) = thisconnection;
        elseif ~isempty(find(funcsynapses(1).ICells==thisconnection(2)));
            funcsynapses.ConnectionsEI(end+1,:) = thisconnection;
        else
            funcsynapses.ConnectionsEUnk(end+1,:) = thisconnection;
        end
    end
end

if ~isempty(funcsynapses.ConnectionsI)
    for a = 1:size(funcsynapses.ConnectionsI)
        thisconnection = funcsynapses.ConnectionsI(a,:);
        if ~isempty(find(funcsynapses.ECells==thisconnection(2)));
            funcsynapses.ConnectionsIE(end+1,:) = thisconnection;
        elseif ~isempty(find(funcsynapses(1).ICells==thisconnection(2)));
            funcsynapses.ConnectionsII(end+1,:) = thisconnection;
        else
            funcsynapses.ConnectionsIUnk(end+1,:) = thisconnection;
        end
    end
end
    

%% If zero lag synapses, categorize unique pairs
if isfield(funcsynapses,'ZeroLag');
    if isfield(funcsynapses.ZeroLag,'CnxnWeightsZ')
        [pre post]=find(funcsynapses.ZeroLag.CnxnWeightsZ>0);
        prepost = unique([pre post],'rows');
        prepost = EliminatePalindromeRows(prepost);
        funcsynapses.ZeroLag.EPairs = prepost;
        
        [pre post]=find(funcsynapses.ZeroLag.CnxnWeightsZ<0);
        prepost = unique([pre post],'rows');
        prepost = EliminatePalindromeRows(prepost);
        funcsynapses.ZeroLag.IPairs = prepost;

        funcsynapses.ZeroLag.ERanges = [];
        funcsynapses.ZeroLag.IRanges = [];        

        epairs = funcsynapses.ZeroLag.EPairs;
        for a  = 1:size(epairs,1)
            if a==1;
                funcsynapses.ZeroLag.ERanges = [];
            end
            temp = CalcZeroLagRange(funcsynapses,epairs(a,1),epairs(a,2),'above');
            if ~isempty(temp)
                funcsynapses.ZeroLag.ERanges(a,:) = temp;
            end
        end
        
        ipairs = funcsynapses.ZeroLag.IPairs;
        for a  = 1:size(ipairs,1)
            if a==1;
                funcsynapses.ZeroLag.ERanges = [];
            end
            temp = CalcZeroLagRange(funcsynapses,ipairs(a,1),ipairs(a,2),'below');
            if ~isempty(temp)
                funcsynapses.ZeroLag.IRanges(a,:) = temp;
            end
        end    
        
    end    
end

if isfield(funcsynapses,'wide')
   %% do same as above basically?? 
end
