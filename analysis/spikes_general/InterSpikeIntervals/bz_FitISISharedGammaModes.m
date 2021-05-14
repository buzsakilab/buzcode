function [GammaFit] = bz_FitISISharedGammaModes(logISIhist,logtimebins,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%   INPUTS
%       logbins         [numtimebins x 1] (or ISIs)
%       logISIhist      [numtimebins x numcells]  probability density (N/(sum*dbin))
%
%   Options
%       'logbase'       (default: 10)
%       'numAS'         number of activated states
%       'figfolder'     a folder to save the figure in
%       'basePath'
%       'AScost_lambda'
%       'AScost_p'
%       'ASguess'       A L1/2 norm cost on Activated states that tries to
%                       push weights to 0 and avoid overlapping AS modes
%
%       'numpad'        (number of bins below/above to pad)
%       'maxNmodes'  
%       'returnNmodes'  'auto' (selects based on drop in error)
%       'Nestimatemethod'  'ascending' or 'descending'
%       'autoNmodes'    'LargeInflection','TSEthresh'
%       'promthresh'    for auto selection (default: 0.01))
%       'showfig'
%       'lambdabounds'
%       'ISIs'          for measuing AIC (maybe just make hist from this?)
%
%   OUTPUTS
%       lambdas: the beta parameter (mean ISI is ks/lambdas)
%       ks: the alpha parameter     (CV is 1/ks)
%% Input Parser

% parse args
p = inputParser;
addParameter(p,'logbase',10)
addParameter(p,'numAS',3)
addParameter(p,'showfig',true)
addParameter(p,'figfolder',false)
addParameter(p,'figname',[])
addParameter(p,'basePath',pwd,@isstr)
addParameter(p,'AScost_lambda',0)
addParameter(p,'AScost_p',2/3)
addParameter(p,'MScost',0)
addParameter(p,'ASguess',[])

addParameter(p,'savecellinfo',false,@islogical)

addParameter(p,'returnNmodes',6)
addParameter(p,'autoNmodes',true)
addParameter(p,'Nestimatemethod','descending')
addParameter(p,'maxNmodes',10)
%addParameter(p,'lambdabounds',[-5 8])
addParameter(p,'logratebounds',[-3.25 3.25])
addParameter(p,'numpad',15)
addParameter(p,'minISIs',300)
addParameter(p,'promthresh',0.05)
addParameter(p,'sequentialreduce',false)
%addParameter(p,'lasso',0)


parse(p,varargin{:})
logbase = p.Results.logbase;
numAS = p.Results.numAS;
SHOWFIG = p.Results.showfig;
figfolder = p.Results.figfolder;
figname = p.Results.figname;
basePath = p.Results.basePath;
AScost_lambda = p.Results.AScost_lambda;
AScost_p = p.Results.AScost_p;
ASguess = p.Results.ASguess;
MScost = p.Results.MScost;

educatedGuess.logrates =  [2.0 2.25 0.9  1.5 1.0 0.5]; %HiGamma Burst Theta LowGamma Irregular10Hz ThetaDouble
educatedGuess.CVs =       [0.2 0.05 0.05 0.3 0.6 0.05];

numpad = p.Results.numpad;
maxNmodes = p.Results.maxNmodes;
returnNmodes = p.Results.returnNmodes;
logratebounds = p.Results.logratebounds; %units: loglambda, e
minISIs = p.Results.minISIs;
sequentialreduce = p.Results.sequentialreduce;
promthresh = p.Results.promthresh;
autoNmodes = p.Results.autoNmodes;
Nestimatemethod = p.Results.Nestimatemethod;
%lasso = p.Results.lasso; %lasso doesn't work because weights sum to 1...
%%

if strcmp(figfolder,'detectionfigures')
    figfolder = [basePath,filesep,'DetectionFigures'];
end
% if length(ISIs)<minISIs
%     return
% end
%% DEV


baseName = bz_BasenameFromBasepath(basePath);

KSfit = false;
if strcmp(logtimebins,'KSfit')
    KSfit = true;
    logtimebins = linspace(-3.5,2.5,100);
    logISIs = log(ISIs);
    logISIhist = hist(logISIs,logtimebins);
end

taubins = logtimebins./log10(exp(1));
logISIhist = logISIhist.* mode(diff(logtimebins))./mode(diff(taubins)); %convert to dtau

sub1msbins = logtimebins<=-2.7;

%zerospkcells = sum(logISIhist,1)==0;

%meanISI = mean(taubins.*logISIhist)./sum(logISIhist);
numcells = size(logISIhist,2);
%% Dev - fake distribution

% numcells = 10;
% GSASstruct.GSlogrates = linspace(-1,0,numcells);
% GSASstruct.GSCVs = ones(1,numcells);
% GSASstruct.GSweights = 0.5.*ones(1,numcells);
% 
% numAS = 2;
% GSASstruct.ASlogrates = [1 2];
% GSASstruct.ASCVs = [0.5 0.1];
% GSASstruct.ASweights  = rand(numcells,numAS);
% 
% [allISIdist1] = GSASmodel((GSASstruct),taubins,numcells,numAS);
% [allISIdist2] = GSASmodel(convertGSASparms(GSASstruct),taubins,numcells,numAS);
% [allISIdist3] = GSASmodel(convertGSASparms(convertGSASparms(GSASstruct),numcells,numAS),taubins,numcells,numAS);
% 
% figure
% subplot(2,2,1)
% imagesc(taubins,[1 numcells],allISIdist1')
% subplot(2,2,2)
% imagesc(taubins,[1 numcells],allISIdist2')
% subplot(2,2,3)
% imagesc(taubins,[1 numcells],allISIdist3')

%% Set up everything for fitting


%Initial Conditions
%init_struct.GSlogrates = -log10(meanISI)-0.5;
init_struct.GSlogrates = -0.5.*ones(1,numcells);
init_struct.GSCVs = 1.5.*ones(1,numcells);
init_struct.GSweights = 0.5.*ones(1,numcells);

if ASguess
    init_struct.ASlogrates = educatedGuess.logrates(1:numAS);
    init_struct.ASCVs = educatedGuess.CVs(1:numAS);
else
    init_struct.ASlogrates = linspace(1,2.5,numAS);
    init_struct.ASCVs = 0.3.*ones(1,numAS);
end
init_struct.ASweights  = 0.5.*ones(numcells,numAS)./(numAS);
init = convertGSASparms(init_struct);

%Upper/Lower Bounds
clear lb ub
lb.GSlogrates = -2.*ones(1,numcells);
lb.GSCVs =      zeros(1,numcells);
lb.GSweights =  zeros(1,numcells);
lb.ASlogrates = 0.3.*ones(1,numAS);
lb.ASCVs =      zeros(1,numAS);
lb.ASweights  = zeros(numcells,numAS);
lb = convertGSASparms(lb);

ub.GSlogrates = 2.*ones(1,numcells);
ub.GSCVs =      4.*ones(1,numcells);
ub.GSweights =  ones(1,numcells);
ub.ASlogrates = 3.*ones(1,numAS);
ub.ASCVs =      2.*ones(1,numAS);
ub.ASweights  = ones(numcells,numAS);
ub = convertGSASparms(ub);

%Make the constraint matrix for all weights to add to 1
Aeq = zeros(numcells,length(ub));
Aeq_ASonly = zeros(numcells,length(ub));
Beq = ones(numcells,1);
for cc = 1:numcells
    thiscell.GSlogrates = zeros(1,numcells);
    thiscell.GSCVs =      zeros(1,numcells);
    thiscell.GSweights =  zeros(1,numcells);
    thiscell.ASlogrates = zeros(1,numAS);
    thiscell.ASCVs =      zeros(1,numAS);
    thiscell.ASweights  = zeros(numcells,numAS);
    thiscell.ASweights(cc,:) = 1;
    Aeq_ASonly(cc,:) = convertGSASparms(thiscell);
    thiscell.GSweights(cc) = 1;
    Aeq(cc,:) = convertGSASparms(thiscell);
end
Aeq_ASonly(Aeq_ASonly~=1)=0;
Aeq(Aeq~=1)=0;

options = optimoptions('fmincon','Algorithm','sqp' ,'UseParallel',false,'Display','iter');%
%try also: 'Algorithm','interior-point''active-set'
%Decrease tolerance.....
options.MaxFunctionEvaluations = 1e8;
options.MaxIterations = 1000; 

%% Fit all the distributions together

if KSfit
costfun = @(GSASparm_vect) ksstat(logISIs,GSASmodel(GSASparm_vect,'sample',numcells,numAS)) ...
    + AScost_lambda.*sum((abs(Aeq_ASonly*GSASparm_vect)).^(AScost_p))...; %L1/2 norm on AS weights to promote sparseness
    + MScost.*sum(sum((logISIhist(sub1msbins,:)-GSASmodel(GSASparm_vect,taubins(sub1msbins),numcells,numAS)).^2)); 
else                                  
costfun = @(GSASparm_vect) sum(sum((logISIhist-GSASmodel(GSASparm_vect,taubins,numcells,numAS)).^2)) ...
    + AScost_lambda.*sum((abs(Aeq_ASonly*GSASparm_vect)).^(AScost_p))...; %L1/2 norm on AS weights to promote sparseness
    + MScost.*sum(sum((logISIhist(sub1msbins,:)-GSASmodel(GSASparm_vect,taubins(sub1msbins),numcells,numAS)).^2)); 
end
fitparms = fmincon(costfun,init,[],[],Aeq,Beq,lb,ub,[],options);
sharedfit = convertGSASparms(fitparms,numcells,numAS);

%% Fit Each cell distribution, starting from the group dists
display('Group fit Complete! Now Fitting each cell independently')
for cc = 1:numcells
    bz_Counter(cc,numcells,'Cell')
    %Initial Conditions
    cinit_struct.GSlogrates = sharedfit.GSlogrates(cc);
    cinit_struct.GSCVs = sharedfit.GSCVs(cc);
    cinit_struct.GSweights = sharedfit.GSweights(cc);

    cinit_struct.ASlogrates = sharedfit.ASlogrates;
    cinit_struct.ASCVs = sharedfit.ASCVs;
    cinit_struct.ASweights  = sharedfit.ASweights(cc,:);
    cinit = convertGSASparms(cinit_struct);

    %Upper/Lower Bounds
    clear clb cub
    clb.GSlogrates = -2.*ones(1,1);
    clb.GSCVs =      zeros(1,1);
    clb.GSweights =  zeros(1,1);
    clb.ASlogrates = 0.*ones(1,numAS);
    clb.ASCVs =      zeros(1,numAS);
    clb.ASweights  = zeros(1,numAS);
    clb = convertGSASparms(clb);

    cub.GSlogrates = 2.*ones(1,1);
    cub.GSCVs =      4.*ones(1,1);
    cub.GSweights =  ones(1,1);
    cub.ASlogrates = 2.5.*ones(1,numAS);
    cub.ASCVs =      2.*ones(1,numAS);
    cub.ASweights  = ones(1,numAS);
    cub = convertGSASparms(cub);

    %Make the constraint matrix for all weights to add to 1
    cAeq = zeros(1,length(cub));
    cAeq_ASonly = zeros(1,length(cub));
    cBeq = ones(1,1);
    thiscell.GSlogrates = zeros(1,1);
    thiscell.GSCVs =      zeros(1,1);
    thiscell.GSweights =  zeros(1,1);
    thiscell.ASlogrates = zeros(1,numAS);
    thiscell.ASCVs =      zeros(1,numAS);
    thiscell.ASweights  = ones(1,numAS);
    cAeq_ASonly = convertGSASparms(thiscell)';
    thiscell.GSweights =  ones(1,1);
    cAeq = convertGSASparms(thiscell)';
    cAeq_ASonly(cAeq_ASonly~=1)=0;
    cAeq(cAeq~=1)=0;
    
    options = optimoptions('fmincon','Algorithm','sqp','UseParallel',false,'Display','off');%
    %try also: 'Algorithm','active-set', 'sqp'
    %Decrease tolerance.....
    options.MaxFunctionEvaluations = 1e8;
    options.MaxIterations = 1000; 

    %% Fit the single-cell distribution 
    thisdist = logISIhist(:,cc);
                                    
    cdifffun = @(GSASparm_vect) sum(sum((thisdist-GSASmodel(GSASparm_vect,taubins,1,numAS)).^2)) ...
    + AScost_lambda.*sum((abs(cAeq_ASonly*GSASparm_vect)).^(AScost_p))... ; %L2/3 norm on AS weights to promote sparseness
    + MScost.*sum(sum((thisdist(sub1msbins)-GSASmodel(GSASparm_vect,taubins(sub1msbins),1,numAS)).^2)); 

    fitparms_singlecell = fmincon(cdifffun,cinit,[],[],cAeq,cBeq,clb,cub,[],options);
%     if(zerospkcells(cc))
%         fitparms_singlecell = nan(size(fitparms_singlecell));
%     end
    GammaFit.singlecell(cc) = convertGSASparms(fitparms_singlecell,1,numAS);
    
    
end
%GammaFit.singlecell(zerospkcells) = [];
GammaFit.taubins = taubins;
GammaFit.ISIdists = logISIhist;
GammaFit.sharedfit = sharedfit;
GammaFit.logtimebins = logtimebins;
GammaFit.numcells = numcells;
GammaFit.parms = p.Results;
GammaFit.initialconditions = init_struct;
%Collapse the structure
singlecell_all = bz_CollapseStruct(GammaFit.singlecell,1);

%%
if SHOWFIG | figfolder
    
GScolor = [0.6 0.4 0];

figure
hist(singlecell_all.GSweights)
%%
[~,sortGSrate] = sort(sharedfit.GSlogrates);
for aa = 1:numAS
    [~,sortASweight{aa}] = sort(sharedfit.ASweights(:,aa));
end

%% ISI dist
meanISIdist = mean(logISIhist,2);
%%
%sortrate = sort(ISIStats.summstats.WAKEstate.meanrate);
fitISI = GSASmodel(sharedfit,taubins,numcells,numAS);
figure
subplot(2,2,1)
imagesc(logtimebins,[1 numcells],logISIhist(:,sortGSrate)')
hold on
plot(logtimebins,-bz_NormToRange(meanISIdist,0.3)+numcells,'k','linewidth',2)
%colorbar
title(figname)
xlim([-3 2])

subplot(2,2,2)
imagesc(logtimebins,[1 numcells],fitISI(:,sortGSrate)')
%colorbar
xlim([-3 2])

subplot(2,2,2)



subplot(2,2,3)
plot(-sharedfit.ASlogrates,log10(sharedfit.ASCVs),'o')
hold on
% scatter(-sharedfit.GSlogrates,sharedfit.GSCVs,...
%     5.*singlecell_all.GSweights+0.00001,log10(ISIStats.summstats.WAKEstate.meanrate(ISIStats.sorts.WAKEstate.ratepE)),...
%     'filled')
scatter(-sharedfit.GSlogrates,log10(sharedfit.GSCVs),...
    10.*singlecell_all.GSweights+0.00001,...
    'filled')
for aa = 1:numAS
scatter(-singlecell_all.ASlogrates(:,aa),log10(singlecell_all.ASCVs(:,aa)),25.*singlecell_all.ASweights(:,aa)+0.00001,...
    'filled')
end
%axis tight
box off
plot(logtimebins([1 end]),[0 0],'k--')
xlabel('Mean ISI (s)');ylabel('CV')
xlim([-3 2])
LogScale('x',10)

subplot(2,2,4)
plot(singlecell_all.GSlogrates,sharedfit.GSlogrates,'.')
xlabel('Single-cell GS rate');ylabel('Group Fit GS rate')

if figfolder
    NiceSave(['GammaModes',figname],figfolder,baseName);
end

%%
figure
subplot(3,3,1)
imagesc(taubins,[1 numcells],logISIhist(:,sortGSrate)')
%colorbar
for aa = 1:numAS
subplot(3,3,3+aa)
imagesc(taubins,[1 numcells],logISIhist(:,sortASweight{aa})')
end


%% Single cell example(s)
numex = 3;
excells = randi(numcells,numex);
figure
for ee = 1:numex
    excell = excells(ee);
subplot(3,3,ee)
    plot(logtimebins,logISIhist(:,excell),'color',[0.5 0.5 0.5],'linewidth',2)
    hold on
    plot(logtimebins,fitISI(:,excell),'k--','linewidth',1)
    plot(GammaFit.logtimebins,...
        GSASmodel(GammaFit.singlecell(excell),...
        GammaFit.taubins),...
        'k','linewidth',2)
    hold on
    plot(logtimebins,...
        LogGamma(singlecell_all.GSlogrates(excell),...
        singlecell_all.GSCVs(excell),...
        singlecell_all.GSweights(excell)',taubins'),'color',GScolor);
    for aa = 1:numAS
        plot(logtimebins,...
            LogGamma(singlecell_all.ASlogrates(excell,aa),...
            singlecell_all.ASCVs(excell,aa),...
            singlecell_all.ASweights(excell,aa)',taubins'),'k');
    end
    box off
    axis tight
    %title(

subplot(3,3,3+ee)
scatter(-singlecell_all.ASlogrates(excell,:),log10(singlecell_all.ASCVs(excell,:)),...
    100*singlecell_all.ASweights(excell,:)+0.00001,'k','filled')
hold on
scatter(-singlecell_all.GSlogrates(excell),log10(singlecell_all.GSCVs(excell)),...
    100*singlecell_all.GSweights(excell)+0.00001,GScolor,'filled')
plot(logtimebins([1 end]),[0 0],'k--')
ylabel('CV');xlabel('mean ISI (s)')
xlim(logtimebins([1 end]));ylim([-2 1])
LogScale('x',10,'exp',true)

end
if figfolder
    NiceSave(['CellExample',figname],figfolder,baseName);
end

end

%%
   
%    
% meanISI = ks./lambdas;
% %Initialize parms
% init = [linspace(-1.5,5.5,returnNmodes)';...    %Lambda 
%     -0.3.*ones(returnNmodes,1);         %K  (used to be CV=0.8...)
%     ones(returnNmodes,1)./(returnNmodes)];             %Weights (normalize later)
%     
% figure
% 
% subplot(4,2,1)
% plot(timebins,logISIhist,'color',[0.5 0.5 0.5],'linewidth',2)
% hold on
% plot(timebins,multigamfun(fitparms{returnNmodes},taubins),'k','linewidth',2)
% for mm = 1:returnNmodes
%     plot(timebins,multigamfun(fitparms{returnNmodes}(returnNmodes.*[0;1;2]+mm),taubins),'r')
% end
% axis tight
% box off
% LogScale('x',logbase)
% 
% subplot(4,2,3)
% stem(log10(1./meanISI),weights)
% LogScale('x',10)
% xlabel('Rate (Hz)')
% ylabel('weight')
% 
% subplot(2,2,2)
% scatter(log10(1./meanISI),1./ks,10)
% %crameri('bilbao')
% LogScale('x',10)
% xlabel('Rate (Hz)');ylabel('1/k (CV)')
% 
% subplot(6,3,11)
% plot(trymodes,errordrop,'o-')
% hold on
% plot(trymodes(putNmodes),errordrop(putNmodes),'r^')
% xlabel('Number of Modes')
% ylabel('dTSE')
% 
% subplot(6,3,12)
% plot(trymodes(1:end-1),inflection,'o-')
% hold on
% plot([0 trymodes(end)],[0 0],'r--')
% plot(trymodes(putNmodes),errordrop(putNmodes),'r^')
% xlabel('Number of Modes')
% ylabel('inflect.')
% %LogScale('y',10)
% 
% subplot(6,3,10)
% plot(trymodes,log10(fiterror),'o-')
% xlabel('Number of Modes')
% ylabel('logTSE')
% %LogScale('y',10)
% 
% if ~sequentialreduce
% subplot(4,2,3)
% plot(timebins,logISIhist,'k','linewidth',2)
% hold on
% plot(timebins,multigamfun(init,taubins),'r','linewidth',2)
% for mm = 1:returnNmodes
%     plot(timebins,multigamfun(init(returnNmodes.*[0;1;2]+mm),taubins),'r')
% end
% title('Initalization')
% end
% 
% for nn = 1:6
% subplot(6,3,nn+12)
% plot(timebins,logISIhist,'k','linewidth',2)
% hold on
% plot(timebins,multigamfun(ICs{nn},taubins),'r--','linewidth',1)
% plot(timebins,multigamfun(fitparms{nn},taubins),'r','linewidth',2)
% for mm = 1:nn
%     plot(timebins,multigamfun(fitparms{nn}(nn.*[0;1;2]+mm),taubins),'r')
% end
% 
% axis tight
% box off
% LogScale('x',logbase)
% 
% end
% 
% end

end
