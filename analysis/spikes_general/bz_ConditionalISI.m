function [ConditionalISI] = bz_ConditionalISI(spikes,conditionalvariable,varargin)
%[ConditionalISI] = ConditionalISI(spikes,conditionalvariable,<options>)
%    
%   INPUTS
%       spikes  cell array (from spikes.times)
%       conditionalvariable struct with .data and .timestamps
%
%Future: multiple conditional variables... (fit together)
%%
p = inputParser;
addParameter(p,'ints',[])
addParameter(p,'normtype','percentile')
addParameter(p,'minX',50)
addParameter(p,'numISIbins',120)
addParameter(p,'ISIbounds',[0.001 100])
addParameter(p,'numXbins',10)
addParameter(p,'Xwin',[0 1])
addParameter(p,'GammaFitParms',[])
addParameter(p,'GammaFit',false)
addParameter(p,'MutInf',true)
addParameter(p,'ISIDist',true)
addParameter(p,'showfig',true)
addParameter(p,'figname',[])
addParameter(p,'basePath',pwd,@isstr)
addParameter(p,'figfolder',false)
addParameter(p,'Xbinoverlap',1)
addParameter(p,'spikecondition',true)


parse(p,varargin{:})
ints = p.Results.ints;
normtype = p.Results.normtype;
minX = p.Results.minX;
numISIbins = p.Results.numISIbins;
logISIbounds = log10(p.Results.ISIbounds);
numXbins = p.Results.numXbins;
Xwin = p.Results.Xwin;
GFParms = p.Results.GammaFitParms;
DO_GammaFit = p.Results.GammaFit;
SHOWFIG = p.Results.showfig;
figname = p.Results.figname;
basePath = p.Results.basePath;
figfolder = p.Results.figfolder;
DO_ISIDist = p.Results.ISIDist;
DO_MutInf = p.Results.MutInf;
Xbinoverlap = p.Results.Xbinoverlap;
spikecondition = p.Results.spikecondition;
%%
baseName = bz_BasenameFromBasepath(basePath);

%%
if iscell(spikes)
    [all_ConditionalISI] = cellfun(@(spk) bz_ConditionalISI(spk,conditionalvariable,varargin{:}),...
        spikes,'UniformOutput',false);
    
    all_ConditionalISI = cat(1,all_ConditionalISI{:});
    ConditionalISI = bz_CollapseStruct( all_ConditionalISI,3,'justcat',true);
    
    %HERE! Fit all cells simultaneously? (ugh, but with separate
    %parameters...)
    return
end

%% Calculate ISIs
ISIs.n = diff(spikes);
ISIs.times = spikes(2:end-1);
ISIs.np1 = ISIs.n(2:end);
ISIs.n = ISIs.n(1:end-1);

%% If Conditional Variable is Intervals - make a state vector
if strcmp(ints,'input')
    %Structure with states.NAMEstate format
    [ conditionalvariable ] = bz_INTtoIDX(conditionalvariable,'sf',100);
    conditionalvariable.data = conditionalvariable.states;
    
    %One bin per state
    ConditionalISI.states = conditionalvariable.statenames;
    numXbins = length(ConditionalISI.states);
    Xwin = [0.5 numXbins+0.5];
    
    
    %Convert to all intervals to restrict spikes...
    conditionalvariable.states = conditionalvariable.states>0;
    conditionalvariable.statenames = {'AllInts'};
    ints = bz_IDXtoINT(conditionalvariable);
    ints = ints.AllIntsstate;
end
%% Restrict ISIs and conditional variable to intervals
conditionalvariable.dt = mode(diff(conditionalvariable.timestamps));

if ~isempty(ints)
    inints = InIntervals(conditionalvariable.timestamps,ints);
    conditionalvariable.data = conditionalvariable.data(inints);
    conditionalvariable.timestamps = conditionalvariable.timestamps(inints);

    inints = InIntervals(ISIs.times,ints);
    ISIs.n = ISIs.n(inints);
    ISIs.np1 = ISIs.np1(inints);
    ISIs.times = ISIs.times(inints);
end

if ~strcmp(normtype,'none')
    [conditionalvariable.data] = NormToInt(conditionalvariable.data,normtype);
end


ISIs.x = interp1(conditionalvariable.timestamps,conditionalvariable.data,ISIs.times,'nearest');
%% Conditional Distribution/Rate, and MutInfo
if DO_ISIDist
    if DO_GammaFit || spikecondition
        [ ConditionalISI.Dist ] = ConditionalHist( [ISIs.x;ISIs.x],log10([ISIs.n;ISIs.np1]),...
            'Xbounds',Xwin,'numXbins',numXbins,'Ybounds',logISIbounds,'numYbins',numISIbins,'minX',minX,...
            'Xbinoverlap',Xbinoverlap);

        %Convert to prob density
        ConditionalISI.Dist.pYX = ConditionalISI.Dist.pYX./mode(diff(ConditionalISI.Dist.Xbins));
    else
        [ ConditionalISI.Dist ] = ConditionalHist( [ISIs.x;ISIs.x],log10([ISIs.n;ISIs.np1]),...
            'Xbounds',Xwin,'numXbins',numXbins,'Ybounds',logISIbounds,'numYbins',numISIbins,'minX',minX,...
            'Xbinoverlap',Xbinoverlap,'conditionby',conditionalvariable.data);
    end
    
    ConditionalISI.Dist.Xocc = hist(conditionalvariable.data,ConditionalISI.Dist.Xbins);
    ConditionalISI.Dist.Xocc = movsum(ConditionalISI.Dist.Xocc,Xbinoverlap.*2-1);
    ConditionalISI.Dist.SpikeRate = ConditionalISI.Dist.Xhist./(ConditionalISI.Dist.Xocc.*conditionalvariable.dt.*2);


end

if DO_MutInf
    ConditionalISI.MutInf = mutualinfo(log10([ISIs.n;ISIs.np1]),[ISIs.x;ISIs.x]);
end

if DO_GammaFit
    %% Set up everything for fitting


    %Initial Conditions
    %init_struct.GSlogrates = -log10(meanISI)-0.5;
    init_struct.GSlogrates = GFParms.GSlogrates.*ones(1,numXbins);
    init_struct.GSCVs = GFParms.GSCVs.*ones(1,numXbins);
    init_struct.GSweights = GFParms.GSweights.*ones(1,numXbins);

    % if ASguess
    %     init_struct.ASlogrates = educatedGuess.logrates(1:numAS);
    %     init_struct.ASCVs = educatedGuess.CVs(1:numAS);
    % else
        init_struct.ASlogrates = GFParms.ASlogrates;
        init_struct.ASCVs = GFParms.ASCVs;
    % end
    init_struct.ASweights  = repmat(GFParms.ASweights,numXbins,1);
    init = convertGSASparms(init_struct);



    %%

    taubins = ConditionalISI.Dist.Ybins./log10(exp(1));
    logISIhist = ConditionalISI.Dist.pYX'.* mode(diff(ConditionalISI.Dist.Xbins))./mode(diff(taubins)); %convert to dtau
    numXbins = numXbins;
    numAS = length(GFParms.ASweights);

    %%


    %Upper/Lower Bounds
    clear lb ub
    lb.GSlogrates = -2.*ones(1,numXbins);
    lb.GSCVs =      zeros(1,numXbins);
    lb.GSweights =  zeros(1,numXbins);
    lb.ASlogrates = 0.3.*ones(1,numAS);
    lb.ASCVs =      zeros(1,numAS);
    lb.ASweights  = zeros(numXbins,numAS);
    lb = convertGSASparms(lb);

    ub.GSlogrates = 2.*ones(1,numXbins);
    ub.GSCVs =      4.*ones(1,numXbins);
    ub.GSweights =  ones(1,numXbins);
    ub.ASlogrates = 3.*ones(1,numAS);
    ub.ASCVs =      2.*ones(1,numAS);
    ub.ASweights  = ones(numXbins,numAS);
    ub = convertGSASparms(ub);

    %Make the constraint matrix for all weights to add to 1
    Aeq = zeros(numXbins,length(ub));
    Aeq_ASonly = zeros(numXbins,length(ub));
    Beq = ones(numXbins,1);
    for cc = 1:numXbins
        thiscell.GSlogrates = zeros(1,numXbins);
        thiscell.GSCVs =      zeros(1,numXbins);
        thiscell.GSweights =  zeros(1,numXbins);
        thiscell.ASlogrates = zeros(1,numAS);
        thiscell.ASCVs =      zeros(1,numAS);
        thiscell.ASweights  = zeros(numXbins,numAS);
        thiscell.ASweights(cc,:) = 1;
        Aeq_ASonly(cc,:) = convertGSASparms(thiscell);
        thiscell.GSweights(cc) = 1;
        Aeq(cc,:) = convertGSASparms(thiscell);
    end
    Aeq_ASonly(Aeq_ASonly~=1)=0;
    Aeq(Aeq~=1)=0;

    options = optimoptions('fmincon','Algorithm','sqp' ,'UseParallel',false,'Display','none');%
    %try also: 'Algorithm','interior-point''active-set'
    %Decrease tolerance.....
    options.MaxFunctionEvaluations = 1e8;
    options.MaxIterations = 1000; 

    %% Fit all the distributions together
    AScost_lambda = 0;
    AScost_p = 1/2;
    MScost = 3;
    sub1msbins = ConditionalISI.Dist.Ybins<=-2.7;

    costfun = @(GSASparm_vect) sum(sum((logISIhist-GSASmodel(GSASparm_vect,taubins,numXbins,numAS)).^2)) ...
        + AScost_lambda.*sum((abs(Aeq_ASonly*GSASparm_vect)).^(AScost_p))...; %L1/2 norm on AS weights to promote sparseness
        + MScost.*sum(sum((logISIhist(sub1msbins,:)-GSASmodel(GSASparm_vect,taubins(sub1msbins),numXbins,numAS)).^2)); 

    fitparms = fmincon(costfun,init,[],[],Aeq,Beq,lb,ub,[],options);
    ConditionalISI.GammaModes = convertGSASparms(fitparms,numXbins,numAS);

    %% Mode Correlations
    [ConditionalISI.GammaModes.GSCorr,ConditionalISI.GammaModes.GScorr_p] = corr(ConditionalISI.Dist.Xbins',ConditionalISI.GammaModes.GSweights','type','Pearson');
    [ConditionalISI.GammaModes.ASCorr,ConditionalISI.GammaModes.AScorr_p] = corr(ConditionalISI.Dist.Xbins',ConditionalISI.GammaModes.ASweights,'type','Pearson');


    %%
    ConditionalISI.GammaModes.GS_R = [ones(size(ConditionalISI.Dist.Xbins')) ConditionalISI.Dist.Xbins']\ConditionalISI.GammaModes.GSweights';
    ConditionalISI.GammaModes.GS_R = ConditionalISI.GammaModes.GS_R(2);
    ConditionalISI.GammaModes.AS_R = [ones(size(ConditionalISI.Dist.Xbins')) ConditionalISI.Dist.Xbins']\ConditionalISI.GammaModes.ASweights;
    ConditionalISI.GammaModes.AS_R = ConditionalISI.GammaModes.AS_R(2,:);
end
%%

if SHOWFIG
    
    
figure
if DO_GammaFit
    GScolor = [0.6 0.4 0];

        testmodel = GSASmodel(ConditionalISI.GammaModes,taubins,numXbins,numAS);
        subplot(2,2,1)
        imagesc(ConditionalISI.Dist.Xbins,ConditionalISI.Dist.Ybins,testmodel)
        hold on
        plot(ConditionalISI.Dist.Xbins,-ConditionalISI.GammaModes.GSlogrates,'ro')
        plot(ConditionalISI.Dist.Xbins([1 end]),-GFParms.GSlogrates.*[1 1],'r--')
        %colorbar
        caxis([0 0.4])
        
        subplot(2,2,3)
        %plot(sharedfit.ASweights
        plot(ConditionalISI.GammaModes.ASlogrates(ConditionalISI.GammaModes.AScorr_p<=0.05),...
            ConditionalISI.GammaModes.AS_R(ConditionalISI.GammaModes.AScorr_p<=0.05),'.k','markersize',10)
        hold on
        plot(ConditionalISI.GammaModes.ASlogrates(ConditionalISI.GammaModes.AScorr_p>0.05),...
            ConditionalISI.GammaModes.AS_R(ConditionalISI.GammaModes.AScorr_p>0.05),'.k','markersize',10)
        %hold on
        plot(GFParms.GSlogrates,ConditionalISI.GammaModes.GS_R,'.','markersize',20,'color',GScolor)
        plot(xlim(gca),[0 0],'k--')
        LogScale('x',10)
        xlabel('Rate (Hz)')
        ylabel('Weight Correlation')
        box off

        subplot(2,2,4)
        plot(ConditionalISI.Dist.Xbins,(ConditionalISI.GammaModes.GSweights),'.','markersize',20,'color',GScolor)
        xlabel('Power');ylabel('P_G_S')
        ylim([0 1])
end


subplot(2,2,2)
yyaxis left
imagesc(ConditionalISI.Dist.Xbins,ConditionalISI.Dist.Ybins,ConditionalISI.Dist.pYX')
bounds = ylim(gca);
hold on
LogScale('y',10)
ylabel('ISI (s)')
bz_AddRightRateAxis
%colorbar




if figfolder
    NiceSave(['ConditionalISI',figname],figfolder,baseName);
end

end
end

