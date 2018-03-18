function [ output_args ] = bz_LFPPowerDist( LFP,varargin )
% bz_LFPPowerDist(LFP) calculates the power distribution of an LFP signal.
%   NOTE: THIS FUNCTION IS UNDER DEVELOPMENT. Feel free to chip in.

%
%INPUTS
%   LFP     [nt x 1] vector of the LFP signal -or- basename
%   (optional)
%   sf          sampling frequency of the LFP (default 1250Hz)
%   int         restrict analysis to specific time intervals 
%   spectype    'wavelet' or 'FFT'
%   frange      [lowf highf]
%   nfreqs      number of frequencies to look at
%   SHOWFIG     true/false (default: false)
%   figfolder   folder to save output figures
%   
%% DEV
%datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/BWData/';
%figfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/misc/PowerDistribution/';
%recname = '20140526_277um';

%load(fullfile(datasetfolder,recname,[recname,'_LFP.mat']))
%load(fullfile(datasetfolder,recname,[recname,'_SleepScore.mat']))
%%
sf = 1250;
%int = StateIntervals.NREMpacket;

%LFP = LFP.CTX;
%LFP = NormToInt(LFP,int,sf,'modZ');

spectype = 'FFT';

%%
warning('This function is still under development. No promises')


frange = [1 128];
nfreqs = 100;

switch spectype
    case 'FFT'
        freqlist = logspace(log10(frange(1)),log10(frange(2)),nfreqs);
        window = 1;
        noverlap = 0.8;
        window = window*sf;
        noverlap = noverlap*sf;
        [spec,freqs,t_FFT] = spectrogram(LFP,window,noverlap,freqlist,sf);
        spec = abs(spec)';
        
        [~,inintIDX] = RestrictInts(t_FFT',int);
        intspec = spec(inintIDX,:);
    case 'wavelet'
        downsamplefactor = 2;
        sf_down = sf./downsamplefactor;
        intLFP = IsolateEpochs2(downsample(LFP,2),int,0,sf_down);
        ncyc = 5;
        [freqlist,t,intspec] = bz_WaveSpec(intLFP,frange,nfreqs,ncyc,1/sf_down,'log');
        intspec = cat(2,intspec{:});
        intspec = abs(intspec)';
end

    
%%
numpowerbins = 200;

%minpower = min(intspec(:)); maxpower = max(intspec(:));
switch spectype
    case 'FFT'
        minpower = 0.5;maxpower = 2.85e3;
    case 'wavelet'
        minpower = -2;maxpower = 2;
end

powerbins = linspace(log10(minpower),log10(maxpower),numpowerbins);
[powerdist_mean] = hist(log10(intspec),powerbins);

%%
figure
%subplot(2,2,1)
imagesc(log2(freqlist),powerbins,powerdist_mean)
axis xy
LogScale('x',2)

%%
alldists = {'birnbaumsaunders','exponential','extreme value','gamma',...
    'generalized extreme value','generalized pareto','inverse gaussian',...
    'logistic','loglogistic','lognormal','nakagami','normal','rayleigh',...
    'rician','tlocationscale','weibull'};

%Distirbutions were removed that consistently showed bad fit to cortical
%LFP data during NREM
testdists = {'gamma','loglogistic','lognormal',...
    'rayleigh','weibull'};

showexamples = [];

D = {};PF = {};bestfit={};
for ff = 1:nfreqs
    ff
    if ismember(ff,showexamples)
        [distfits] = allfitdist(intspec(:,ff),'PDF');
    else
        [distfits] = allfitdist(intspec(:,ff));
    end
    
    %Keep only the distributions tested
    distnames = {distfits(:).DistName};
    keepdists = ismember(distnames,testdists);
    distfits = distfits(keepdists);
    
    %Check the best and worst-fitting distribution
    bestfit{ff} = distfits(1).DistName;
    worstfit{ff} = distfits(end).DistName;
    
    %Sort alphabetically
    distnames = {distfits(:).DistName};
    [distnames,sortdist] = sort(distnames);
  
    D{ff} = distfits(sortdist);
end

%%
AICs = cellfun(@(X) [X(:).AIC],D,'UniformOutput',false);
AICs = cat(1,AICs{:});

ndists = length(distnames);
bestfitmat = zeros(nfreqs,ndists);
for nn = 1:ndists
    bestfitmat(:,nn) = nn.*strcmp(distnames{nn},bestfit);
end
bestfitmat(bestfitmat==0)=nan;
%%
figure
subplot(2,2,1)
    imagesc(log2(freqlist),powerbins,powerdist_mean)
    axis xy
    LogScale('x',2);LogScale('y',10)
    xlabel('f (Hz)')
    ylabel('Power (AU)')
    title([spectype,' Power Distribution'])
    
subplot(4,2,5)
    plot(log2(freqlist),log10(AICs),'LineWidth',1)

    legend(distnames,'location','southwest')
    axis tight
    LogScale('x',2);%LogScale('y',10)
    xlabel('f (Hz)')
    box off
    ylabel('log(AIC)')
    
subplot(6,2,11)
        plot(log2(freqlist),bestfitmat,'.','markersize',20)
        set(gca,'YTick',1:ndists)
        set(gca,'YTickLabels',distnames)
    LogScale('x',2);%LogScale('y',10)
    xlabel('f (Hz)')
    box off
    
    NiceSave('LFPPowerDist',figfolder,recname)



end

