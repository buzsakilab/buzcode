function ccgjitteroutput = CCG_jitter_all(spikesinput,shankstouse,plotting,customoutputname)
% Looks for connections between all cells in a given set of recordings with
% the input basename.  Calls Shige+Asohan's function ccg_jitter.m (modified
% by me.  Saves output to structure ccgjitteroutput.
% Brendon Watson 2012

%% gathering whether we'll plot, default is no plot, if user did not enter it.
if ~exist('shankstouse','var')
    shankstouse = [];
end

if ~exist('plotting','var');
    plotting = 0;
end

inputtype = whos('spikesinput');
switch inputtype.class
    case 'tsdArray'
        [spiket,spikeind] = oneSeries(spikesinput);
        spiket = TimePoints(spiket);
        samplerate = 10000;%%Based on usual TSD object sample rate<< NEED TO GENERALIZE THIS, HAVEN'T FIGURED OUT HOW YET
        numclus = size(spikesinput,1);
    case 'char'
        [spiket, spikeind, numclus, iEleClu, spikeph] = ReadEl4CCG2(spikesinput,shankstouse);
        tempPar = LoadPar([basename,'.xml']);
        samplerate = tempPar.SampleRate;
        clear tempPar
end

ccgjitteroutput(1,1).SpikeIndices = spikeind;


%% Read in basic cluster and spiking info, also recording info from hard drive

%% Loop through each pair of cells, record output of shige's function into my struct
for a=1:numclus      
    for b = (a+1):numclus
%         for a = 1:15  
%             for b = (a+1):15;
        
        [ccgjitteroutput(a,b).ccgR,... 
        ccgjitteroutput(a,b).tR,...
        ccgjitteroutput(a,b).GSPExc,...
        ccgjitteroutput(a,b).GSPInh,...
        ccgjitteroutput(a,b).ccgjMtx,...
        ccgjitteroutput(a,b).ccgjstats,...
        ] = CCG_jitter(spiket,spikeind,a,b,samplerate,samplerate/1000,30,3,1000,0.01,0);

        disp(['finished pair including ',num2str(a),' & ',num2str(b)])
    end
end

%% Classify E, I connections, cells and zerolag
ccgjitteroutput = ccg_jitter_classify(ccgjitteroutput);

%% plot outputs using other functions,  if the user chose do to so, default is no plot
if plotting
    try 
        CCG_jitterbasename_plotpositives(ccgjitteroutput);
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
if exist('customoutputname','var')
    save(customoutputname,'ccgjitteroutput');
elseif strcmp(inputtype.class,'char')
    save([spikesinput '_ccgjitteroutput.mat'],'ccgjitteroutput');
else
    save(['ccgjitteroutput.mat'],'ccgjitteroutput')
end