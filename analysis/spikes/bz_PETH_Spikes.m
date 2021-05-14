function PETHSpikes = bz_PETH_Spikes(events,varargin)

%% input parsing

p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'cellnums',[],@isvector);
addParameter(p,'eventCategory','timestamps',@isstr);
addParameter(p,'eventNumber',1,@isvector);
addParameter(p,'secondsBefore',0.5,@isvector);
addParameter(p,'secondsAfter',1,@isvector);
addParameter(p,'binWidth',0.1,@isvector);
addParameter(p,'plotting',0.1,@isvector);
addParameter(p,'saveMat',false,@islogical);

parse(p,varargin{:})

basepath = p.Results.basepath;
cellnums = p.Results.cellnums;
eventCategory = p.Results.eventCategory;
eventNumber = p.Results.eventNumber;
secondsBefore = p.Results.secondsBefore;
secondsAfter = p.Results.secondsAfter;
binWidth = p.Results.binWidth;
plotting = p.Results.plotting;
saveMat = p.Results.saveMat;

basename = bz_BasenameFromBasepath(basepath);

%% Events 
% events can have one of two formats: buzcode events (name given) or
% timestamps (in seconds)

if ischar(events) %presumped to specify an xxx.events.mat type file
    eventpath = fullfile(basepath,[basename,'.',events,'.events.mat']);
    if ~exist(eventpath,'file')
        error(['No file at: ' eventpath])
        return
    end
    load(eventpath)
    eval(['evs = ' events ';'])
    eval(['triggerTimes = evs. ' eventCategory '(:,eventNumber);'])
elseif isvector(events)
    triggerTimes = events;
end
numtriggers = length(triggerTimes);
    
%% Get spikes
if ~isempty(cellnums)
    spikes = bz_GetSpikes('basepath',basepath,'UID',cellnums);
else
    spikes = bz_GetSpikes('basepath',basepath);
end
times = spikes.times;
numCells = size(times,2);

%% Gather PETHs
secondsBefore = binWidth*floor(secondsBefore/binWidth); %round secondsBefore so it is equally divisible by the binwidth... so there will be a bin starting at the event time
secondsAfter = binWidth*floor(secondsAfter/binWidth);

% loop through trigger events
relativeBins = -secondsBefore:binWidth:secondsAfter;
numBins = length(relativeBins)-1;

counts = zeros(numtriggers,numBins,numCells);
for evn = 1:numtriggers
    for cen = 1:numCells
        thiscellthiseventcounts = histcounts(times{cen},[triggerTimes(evn)+relativeBins]);
%         thiscellthiseventcounts = reshape(thiscellthiseventcounts,[1,1,numbins]);%rehape to add to 
        counts(evn,:,cen) = counts(evn,:,cen) + thiscellthiseventcounts;
    end
end

%% Prepare final output
parameters = p.Results;
parameters.baseName = basename;
parameters.secondsBefore = secondsBefore;
parameters.secondsAfter = secondsAfter;

dimensions = v2struct(numtriggers,numCells,numBins);

PETHSpikes = v2struct(counts,relativeBins,dimensions,parameters);


%% save data
if saveMat    
    if isstr(events)
        [file,path] = uiputfile('.mat','Specify save file name',[basename,'.PETHSpikesTrigBy',events,'.mat']);
    else
        [file,path] = uiputfile('.mat','Specify save file name',[basename,'.PETHSpikesTrigByEvents.mat']);
    end
    
    save(fullfile(path,file),'PETHSpikes')
end

%% plotting
if plotting;
    figure;
    imagesc(squeeze(sum(counts,1))')
    ylabel('Cell Number')
    xlabel('Seconds from Triggers')
    title('Average repsonse per cell across all triggers')
    tickplaces = get(gca,'XTick');
    secondticks = tickplaces*binWidth-secondsBefore;
    ticktext = cellstr(num2str(secondticks'));
    set(gca,'XTickLabel',ticktext)
    hold on
    yl = ylim;
    zx = find(secondticks==0);%where zero is
    zx = tickplaces(zx);
    if ~isempty(zx)
        plot([zx zx],yl,'k')
    end
end

end