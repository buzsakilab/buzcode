function [corrmat,sites,h] = LFPPairMat(varargin)
%LFPPairMat(lfpfilename,varargin) calculates the pairwise relation (options
%below) between all LFP channels
%
%INPUTS
%   lfpfilename     ?????    
%   (optional)      put text string followed by value
%       'channels'  only calculate relationship metric between specific
%                   channels
%       'ints'      time intervals during which to calculate sitemap
%       'frange'    frequency range to calculate coherence/correlation
%       'metric'    pairwise metric to use
%           'corr'  correlation between raw signal
%           'powercorr' correlation between power of filtered signal
%           'ISPC'  (default) (functional) AKA phase coherence
%           'filtcorr'  correlation between filtered signal
%       'FMA'
%       'numload'   max number of channels to load at a time, decreases
%                   computing time but increases RAM pressure
%       'randint'   saves time by selecting a random subset of time to
%                   calculate pairwise metric within - default... ???
%       'basepath'  location of files, default is undeclared
%       'basename'  name of recordings, default is undeclared
%
%OUTPUTS
%   corrmat         matrix of pairwise relationships
%   sites           order of recording site indices
%           
%TO DO
%   -insert checkpoint to asses how approprite numload is - based on how
%   much RAM LFPss1 is taking?
%   -clean up metric implementation 
%
%% inputParse for Optional Inputs and Defaults
frange = 'gamma';
a = 1;
while a<length(varargin);
    if isa(varargin{1},'char')
        switch varargin{a}
            case 'frange'
                frange = varargin{a+1};
                varargin(a:a+1) = [];
                a = a-1;
            case 'basepath'
                basepath = varargin{a+1};
                varargin(a:a+1) = [];
                a = a-1;
            case 'basename'
                basename = varargin{a+1};
                varargin(a:a+1) = [];
                a = a-1;
        end    
    end
    a = a+1;
end

p = inputParser;

defaultFMA = true;
checkFMA = @(x) islogical(x);

defaultInts = [0 Inf];
checkInts = @(x) size(x,2)==2 && isnumeric(x);

defaultNumload = 67; %approx number of LFP channels to load at a time

defaultChannels = 'all';
checkChannels = @(x) strcmp(x,'all') || isnumeric(x);

defaultMetric = 'ISPC';
validMetrics = {'filtcorr','corr','powercorr','ISPC'};
checkMetric = @(x) any(validatestring(x,validMetrics));

% defaultFrange = 'gamma';
% validFranges = {'delta','theta','spindles','gamma','ripples'};
% checkFrange = @(x) any(validatestring(x,validFranges)) || size(x) == [1,2];


addParameter(p,'ints',defaultInts,checkInts)
addParameter(p,'metric',defaultMetric,checkMetric)
% addParameter(p,'frange',defaultFrange,checkFrange)
addParameter(p,'channels',defaultChannels,checkChannels);
addParameter(p,'numload',defaultNumload,@isnumeric);

parse(p,varargin{:})

%% FMA Toolbox
%if using FMA Toolbox and A current session hasn't been declared, open an
%FMA session
% if FMA && ~exist('DATA','var')
if exist('basepath','var') && exist('basename','var')
    SetCurrentSession('basepath',basepath,'basename',basename);
else
    SetCurrentSession
end
% end
%%
numload = p.Results.numload/2;
numsites = length(p.Results.channels);
numgroups = ceil(numsites/numload);
sites = p.Results.channels;
groupstarts = round(linspace(0,numsites+1,numgroups+1)); 
groups = discretize(1:numsites,groupstarts);
%% Calculate the pairwise matrix
corrmat = ones(numsites);
%%
for ss1 = 1:numgroups
    %%
    tic
    display(['ss1: ',num2str(ss1),' of ',num2str(numgroups)])
    group1ind = find(groups==ss1);
    %loadLFPss1 - load numload at a time?
    LFPss1 = double(GetLFP(sites(group1ind),'intervals',p.Results.ints));
    
    %filter LFPss1 for power
    if ~strcmp(p.Results.metric,'corr')
%         LFPss1 = FilterLFP(LFPss1,'passband',p.Results.frange);
        LFPss1 = FilterLFP(LFPss1,'passband',frange);

        if strcmp(p.Results.metric,'ISPC')
            [LFPss1,~] = Phase(LFPss1);
        elseif strcmp(p.Results.metric,'powercorr')
            [~,LFPss1] = Phase(LFPss1);
        elseif ~strcmp(p.Results.metric,'filtcorr')
            display(['Doing correlation between filtered signal,'...
                'if this isn"t what you want, you have a problem'])
        end
    end
    
    if strcmp(p.Results.metric,'ISPC')
    	corrmat(group1ind,group1ind) = ISPCmat(LFPss1(:,2:end),LFPss1(:,2:end));
    else
        %get correlation between first numload sites
        corrmat(group1ind,group1ind) = corr(LFPss1(:,2:end),LFPss1(:,2:end),'type','spearman');
    end
        toc
    %%
    for ss2 = (ss1+1):numgroups
        %%
        tic
        display(['ss2: ',num2str(ss2),' of ',num2str(numgroups)])
        group2ind = find(groups==ss2);
        %loadLFPss2 - load numload at a time
        LFPss2 = double(GetLFP(p.Results.channels(group2ind),'intervals',p.Results.ints));
        %filter LFPss2 for power.... make this all case/switch
        if ~strcmp(p.Results.metric,'corr')
%             LFPss2 = FilterLFP(LFPss2,'passband',p.Results.frange);
            LFPss2 = FilterLFP(LFPss2,'passband',frange);
            
            if strcmp(p.Results.metric,'ISPC')
                [LFPss2,~] = Phase(LFPss2);
            elseif strcmp(p.Results.metric,'powercorr')
                [~,LFPss2] = Phase(LFPss2);
            elseif ~strcmp(p.Results.metric,'filtcorr')
                display(['Doing correlation between filtered signal,'...
                    'if this isn"t what you want, you have a problem'])
            end
        end
        
        if strcmp(p.Results.metric,'ISPC')
            corrmat(group1ind,group2ind) = ISPCmat(LFPss1(:,2:end),LFPss2(:,2:end));
        else
            %measure correlation or coherence beween this numload and the other
            %numload
            corrmat(group1ind,group2ind) = corr(LFPss1(:,2:end),LFPss2(:,2:end),'type','spearman');
        end
        corrmat(group2ind,group1ind) = corrmat(group1ind,group2ind)'; %symmetry :D
        toc
        %%
    end
 
   %%

end

h = figure;
imagesc(corrmat)
colorbar



end