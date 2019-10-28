function [SpikeLFPCoupling] = bz_GenSpikeLFPCoupling(spikes,LFP,varargin)
% SpikeLFPCoupling = GenSpikeLFPCoupling(spikes,LFP)
%
%INPUT
%   spikes          structure with fields (from bz_getSpikes)
%                           .times
%                       -or-
%                       {Ncells} cell array of spiketimes for each cell
%   LFP                 structure with fields (from bz_GetLFP)
%                           lfp.data
%                           lfp.timestamps
%                           lfp.samplingRate 
%   (optional)
%       'frange'
%       'tbin'
%       'waveparms'...
%       'nfreqs'
%       'int'
%       'ISIpower'      true: calculate mutual information between power
%                       ISI distribution
%       'cellLFPchannel'  The local LFP channel associated with each cell (if
%                       sites are in different electroanatomical regions)
%       'cellsort'        'pca','none','sortf','sortorder','celltype','rate'
%       'controls'        'thinspikes','jitterspikes','shufflespikes'
%       'downsample'
%       'cellclass'     cell array with string label for class of each cell
%       'jittersig'     true/false for jittered spikes significance
%       'showFig'       true/false
%       'spikeLim'      limit number of spikes to look at for each cell 
%                       (randomly omits spikes, default: Inf)
%
%OUTPUT
%        SpikeLFPCoupling
%           .freqs
%           .pop
%               .popcellind
%               .cellpopidx
%           .cell                all results are [cells x freqs x channels]
%               .ratepowercorr      
%               .ratepowersig       
%               .spikephasemag      
%               .spikephaseangle    
%               .spikephasesig
% 
%           .detectorinfo
%               .detectorname
%               .detectiondate
%               .detectionintervals
%               .detectionchanne;
%
%TO DO
%   -Put synch/spikerate in same section
%   -multiple ints - this (along with controls) will require making
%   subfunctions to loop
%   -clean and buzcode
%
%DLevenstein 2016
%% inputParse for Optional Inputs and Defaults
p = inputParser;


checkInt = @(x) size(x,2)==2 && isnumeric(x) || isa(x,'intervalSet');
checkFrange = @(x) isnumeric(x) && length(x(1,:)) == 2 && length(x(:,1)) == 1;

validSorttypes = {'pca','none','sortf','sortorder','celltype','rate'};
checkSorttype = @(x) any(validatestring(x,validSorttypes)) || size(x) == [1,2];

addParameter(p,'int',[0 Inf],checkInt)
addParameter(p,'frange',[1 128],checkFrange)
addParameter(p,'nfreqs',100,@isnumeric)
addParameter(p,'ncyc',5,@isnumeric)
addParameter(p,'synchdt',0.005,@isnumeric)
addParameter(p,'synchwin',0.02,@isnumeric)
addParameter(p,'sorttype','rate',checkSorttype)
addParameter(p,'DOWNSAMPLE',false,@isnumeric)
addParameter(p,'cellclass',0)
addParameter(p,'channel',[])
addParameter(p,'jittersig',false)
addParameter(p,'showFig',true)
addParameter(p,'saveFig',false)
addParameter(p,'saveMat',false)
addParameter(p,'ISIpower',false)
addParameter(p,'spikeLim',Inf)

parse(p,varargin{:})
%Clean up this junk...
int = p.Results.int;
nfreqs = p.Results.nfreqs;
frange = p.Results.frange;
ncyc = p.Results.ncyc;
synchdt = p.Results.synchdt;
jittersig = p.Results.jittersig;
SHOWFIG = p.Results.showFig;
SAVEMAT = p.Results.saveMat;
figfolder = p.Results.saveFig; 
usechannel = p.Results.channel; 
ISIpower = p.Results.ISIpower; 
spikeLim = p.Results.spikeLim;
subpop = p.Results.cellclass;

%% Deal with input types

if isa(spikes,'cell')
    spikes_temp.numcells = length(spikes);
    spikes_temp.times = spikes;
    spikes = spikes_temp;
    clear spikes_temp
end

if isa(int,'intervalSet')
    int = [Start(int,'s'), End(int,'s')];
end


if ~isempty(usechannel)
    usechannel = ismember(LFP.channels,usechannel);
    LFP.data = LFP.data(:,usechannel);
end

%Downsampling
if p.Results.DOWNSAMPLE
    assert((LFP.samplingRate/p.Results.DOWNSAMPLE)>2*max(frange),'downsample factor is too big...')    
    [ LFP ] = bz_DownsampleLFP(LFP,p.Results.DOWNSAMPLE);
end
%% Subpopulations
if isequal(subpop,0)
        numpop = 1;
        popcellind = {1:numcells};
elseif isequal(subpop,'done')
else
        pops = unique(subpop);
        numpop = length(pops);
        for pp = 1:numpop
            popcellind{pp} = find(ismember(subpop,pops{pp}));
        end
end

cellpopidx = zeros(1,spikes.numcells);


%% Calculate spike matrix
spikemat = bz_SpktToSpkmat(spikes,'binsize',p.Results.synchwin,'dt',synchdt);

inint = InIntervals(spikemat.timestamps,int);
spikemat.data = spikemat.data(inint,:);
spikemat.timestamps = spikemat.timestamps(inint);


%% Calculate ISIs (if needed), Apply interval/spike number limits here

%Take only spike times in intervals
spikes.inint = cellfun(@(X) InIntervals(X,int),spikes.times,'UniformOutput',false);

if ISIpower
    %Calculate (pre/postceding ISI) for each spike
    spikes.ISIs_prev = cellfun(@(X) [nan; diff(X)],spikes.times,'UniformOutput',false);
    spikes.ISIs_next = cellfun(@(X) [diff(X); nan],spikes.times,'UniformOutput', false);
end

%Apply spikelimit here
spikes.toomany = cellfun(@(X) (sum(X)-spikeLim).*((sum(X)-spikeLim)>0),spikes.inint);
spikes.numcells = length(spikes.times);
for cc = 1:spikes.numcells
    spikes.inint{cc}(randsample(find(spikes.inint{cc}),spikes.toomany(cc)))=false;
end

%Restrict to up to spikelim spikes in interval
spikes.times = cellfun(@(X,Y) X(Y),spikes.times,spikes.inint,'UniformOutput',false);

%% Processing LFP - filter etc

%HERE: loop channels
for cc = 1:length(LFP.channels)
    chanID = LFP.channels(cc);
    if length(LFP.channels)>1
         fprintf('Channel %d (%d of %d)\n',chanID,cc,length(LFP.channels))
    end

    switch nfreqs  
        case 1
            %Single frequency band - filter/hilbert
            LFP_filt = bz_Filter(LFP,'passband',frange,'order',ncyc,...
                'intervals',int,'filter','fir1','channels',chanID,'fast',false);
            LFP_filt.data = LFP_filt.hilb; 
            freqs = [];
            clear filtered

            %Normalize Power to Mean Power
            LFP_filt.data = LFP_filt.data./mean(abs(LFP_filt.data));

        otherwise
            %Multiple frequencies: Wavelet Transform
            LFP_filt = bz_WaveSpec(LFP,'intervals',int,'showprogress',true,...
                'ncyc',ncyc,'nfreqs',nfreqs,'frange',frange,'chanID',chanID); 
            freqs = LFP_filt.freqs;

            %Normalize power to mean power for each frequency
            LFP_filt.data = bsxfun(@(X,Y) X./Y,LFP_filt.data,nanmean(abs(LFP_filt.data),1));
    end

    %Get Power/Phase at each spike matrix time point and each spike time
    spikemat.filtLFP = interp1(LFP_filt.timestamps,LFP_filt.data,spikemat.timestamps,'nearest');
    %Get complex-valued filtered LFP at each spike time
    if spikes.numcells>50
        disp('Interpolating LFP at each spike... If this is prohibitive (time or RAM), try using ''spikeLim''')
    end
    for nn = 1:spikes.numcells
        bz_Counter(nn,spikes.numcells,'Interpolating Cell')
        spikes.filtLFP{nn} = interp1(LFP_filt.timestamps,LFP_filt.data,spikes.times{nn},'nearest');
    end

    %% Population Synchrony: Phase Coupling and Rate Modulation
    for pp = 1:numpop
        if length(popcellind{pp}) == 0
            popcoupling.(pops{pp}).powercorr = [];
            popcoupling.(pops{pp}).phasemag = [];
            popcoupling.(pops{pp}).phaseangle = [];
            numpop = numpop-1;
            continue
        end
        cellpopidx(popcellind{pp}) = pp;
        numpopcells = length(popcellind{pp});
        popsynch = sum(spikemat.data(:,popcellind{pp})>0,2)./numpopcells;
        popsynch = popsynch./mean(popsynch);

        %Calculate Synchrony-Power Coupling as correlation between synchrony and power
        [popcoupling.(pops{pp}).powercorr(:,cc)] = corr(popsynch,abs(spikemat.filtLFP),'type','spearman','rows','complete');

        %Synchrony-Phase Coupling (magnitude/angle of power-weighted mean resultant vector)
        resultvect = nanmean(abs(spikemat.filtLFP).*bsxfun(@(popmag,ang) popmag.*exp(1i.*ang),...
            popsynch,angle(spikemat.filtLFP)),1);
        popcoupling.(pops{pp}).phasemag(:,cc) = abs(resultvect);
        popcoupling.(pops{pp}).phaseangle(:,cc) = angle(resultvect);
   
    end

    %% Cell Rate-Power Modulation

    %Spike-Power Coupling
    [ratepowercorr(:,:,cc),ratepowersig(:,:,cc)] = corr(spikemat.data,abs(spikemat.filtLFP),'type','spearman','rows','complete');

    %% Cell Spike-Phase Coupling
    
    %Find filtered LFP at the closest LFP timepoint to each spike.

    
    [spikephasemag_cell,spikephaseangle_cell] = cellfun(@(X) spkphase(X),spikes.filtLFP,...
        'UniformOutput',false);
    spikephasemag(:,:,cc) = cat(1,spikephasemag_cell{:});
    spikephaseangle(:,:,cc) = cat(1,spikephaseangle_cell{:});

    if jittersig
        %Jitter for Significane
        numjitt = 100;
        jitterwin = 2/frange(1);

        jitterbuffer = zeros(numcells,nfreqs,numjitt);
        for jj = 1:numjitt
            if mod(jj,10) == 1
                display(['Jitter ',num2str(jj),' of ',num2str(numjitt)])
            end
            jitterspikes = JitterSpiketimes(spikes.times,jitterwin);
            jitterLFP = cellfun(@(X) interp1(LFP_filt.timestamps,LFP_filt.data,X,'nearest'),...
                jitterspikes,'UniformOutput',false);
    
            phmagjitt = cellfun(@(X) spkphase(X),jitterLFP,'UniformOutput',false);
            jitterbuffer(:,:,jj) = cat(1,phmagjitt{:});
        end
        jittermean = mean(jitterbuffer,3);
        jitterstd = std(jitterbuffer,[],3);
        spikephasesig(:,:,cc) = (spikephasemag-jittermean)./jitterstd;
    end
    %% Example Figure : Phase-Coupling Significance
    % cc = 6;
    % figure
    %     hist(squeeze(jitterbuffer(cc,:,:)),10)
    %     hold on
    %     plot(spikephasemag(cc).*[1 1],get(gca,'ylim')./4,'r','LineWidth',2)
    %     plot(spikephasemag(cc),get(gca,'ylim')./4,'ro','LineWidth',2)
    %     xlabel('pMRL')
    %     ylabel('Number of Jitters')
    %     xlim([0 0.2])

    %% Calculate mutual information between ISI distribution and power
    
    if ISIpower
        for nn = 1:spikes.numcells
            bz_Counter(nn,spikes.numcells,'Calculating Power-ISI MI. Cell')
            [totmutXPow(nn,:,cc)] = ISIpowerMutInf(spikes.ISIs_prev{nn}(spikes.inint{nn}),...
                spikes.ISIs_next{nn}(spikes.inint{nn}),spikes.filtLFP{nn});
        end
        
    end
    
    %%
%     randsample(spikes.numcells,1);
%     figure
%     subplot(2,2,1)
%         imagesc(log2(freqs),[1 10],squeeze(mean(totmutXPow,1))')
%         LogScale('x',2)
%     subplot(2,2,2)
%         imagesc(log2(freqs),[1 10],squeeze(totmutXPow(randsample(spikes.numcells,1),:,:))')
%     %caxis([0 0.2])
%     LogScale('x',2)
end

clear LFP_filt

%% Output

SpikeLFPCoupling.freqs = freqs;
SpikeLFPCoupling.pop = popcoupling;
%SpikeLFPCoupling.pop.popcellind = popcellind; %update to match buzcode cell class input
%SpikeLFPCoupling.pop.cellpopidx = cellpopidx;
SpikeLFPCoupling.cell.ratepowercorr = ratepowercorr;
SpikeLFPCoupling.cell.ratepowersig = ratepowersig;
SpikeLFPCoupling.cell.spikephasemag = spikephasemag;
SpikeLFPCoupling.cell.spikephaseangle = spikephaseangle;
if ISIpower
    SpikeLFPCoupling.cell.ISIpowermodulation = totmutXPow;
end
if exist('spikephasesig','var')
    SpikeLFPCoupling.cell.spikephasesig = spikephasesig;
end

SpikeLFPCoupling.detectorinfo.detectorname = 'bz_GenSpikeLFPCoupling';
SpikeLFPCoupling.detectorinfo.detectiondate = datetime('today');
SpikeLFPCoupling.detectorinfo.detectionintervals = int;
SpikeLFPCoupling.detectorinfo.detectionchannel = LFP.channels;

if SAVEMAT
    save(savefile,'SpikeLFPCoupling')
end

%% Figures
if SHOWFIG
    
    %Sorting (and other plot-related things)
    switch p.Results.sorttype
        case 'pca'
            [~,SCORE,~,~,EXP] = pca(cellpower);
            [~,spikepowersort] = sort(SCORE(:,1));

            [~,SCORE,~,~,EXP_phase] = pca(cellphasemag);
            [~,spikephasesort] = sort(SCORE(:,1));
            sortname = 'PC1';
        case 'none'
            spikepowersort = 1:numcells;
        case 'fsort'
            fidx = interp1(freqs,1:nfreqs,sortf,'nearest'); 
            [~,spikepowersort] = sort(cellpower(:,fidx));
            [~,spikephasesort] = sort(cellphasemag(:,fidx));
            sortname = [num2str(sortf) 'Hz Magnitude'];
        case 'rate'
            spkrt = cellfun(@length,spikes.times);
            [~,spikepowersort] = sort(spkrt);
            spikephasesort = spikepowersort;
            sortname = 'Firing Rate';
        otherwise
            spikepowersort = sortidx;
    end
    
    
    if length(LFP.channels)==1
        
    switch nfreqs  
        case 1
    %% Figure: 1 Freq Band
    figure
        subplot(3,2,1)
            hold on
            hist(ratepowercorr)
            for pp = 1:numpop
                plot([popcoupling.(pops{pp}).powercorr,popcoupling.(pops{pp}).powercorr],...
                    get(gca,'ylim'))
            end
            xlabel('Rate-Power Correlation')
            title('Rate-Power Coupling')
            ylabel('# Cells')
        subplot(3,2,[2 4])
            for pp = 1:numpop
                polar(spikephaseangle(popcellind{pp}),(spikephasemag(popcellind{pp})),'o')
                hold on
                polar([0 popcoupling.(pops{pp}).phaseangle],[0 popcoupling.(pops{pp}).phasemag])

            end
            title('Spike-Phase Coupling')
        subplot(3,2,3)
            for pp = 1:numpop
                hold on
                plot(ratepowercorr(popcellind{pp}),(spikephasemag(popcellind{pp})),'o')
                plot(popcoupling.(pops{pp}).powercorr,popcoupling.(pops{pp}).phasemag,'*')
            end
            %LogScale('y',10)
            xlabel('Rate-Power Correlation')
            ylabel('Spike-Phase Coupling Magnitude')
        
        if figfolder
           NiceSave('SpikeLFPCoupling',figfolder,[]) 
        end

        otherwise
    %% Figure: Spectrum
    posnegcolor = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);

    figure
        subplot(3,2,1)
            hold on
            for pp = 1:numpop
                plot(log2(freqs),popcoupling.(pops{pp}).powercorr)
            end
            plot(log2(freqs([1 end])),[0 0],'k--')
            axis tight
            box off
            LogScale('x',2)
            title('Pop. Synchrony - Power Correlation')
            xlabel('f (Hz)');ylabel('rho')

    %     subplot(3,2,3)
    %         hold on
    %             [hAx,hLine1,hLine2] = plotyy(log2(freqs),cat(1,synchcoupling.phasemag),...
    %                 log2(freqs),mod(cat(1,synchcoupling.phaseangle),2*pi))
    %         LogScale('x',2)
    %         hLine1.LineStyle = 'none';
    %         hLine2.LineStyle = 'none';
    %         hLine1.Marker = 'o';
    %         hLine2.Marker = '.';
    %         title('Pop. Synchrony - Phase Coupling')
    %         xlabel('f (Hz)');
    %         ylabel(hAx(1),'Phase Coupling Magnitude')
    %         ylabel(hAx(2),'Phase Coupling Angle')

        subplot(3,2,3)
            hold on
            for pp = 1:numpop
                plot(log2(freqs),popcoupling.(pops{pp}).phasemag)
            end
            axis tight
            box off
            LogScale('x',2)
            title('Pop. Synchrony - Phase Coupling')
            xlabel('f (Hz)');
            ylabel('Phase Coupling Magnitude')
        subplot(3,2,5)
            hold on
            for pp = 1:numpop
                plot(log2(freqs),popcoupling.(pops{pp}).phaseangle,'o')
                plot(log2(freqs),popcoupling.(pops{pp}).phaseangle+2*pi,'o')
            end
            LogScale('x',2)
            title('Pop. Synchrony - Phase Coupling')
            xlabel('f (Hz)');
            axis tight
            box off
            ylim([-pi 3*pi])
            ylabel('Phase Coupling Angle')


        subplot(2,2,2)
            imagesc(log2(freqs),1:spikes.numcells,ratepowercorr(spikepowersort,:))
            xlabel('f (Hz)');ylabel(['Cell - Sorted by ',sortname])
            LogScale('x',2)
            title('Rate - Power Correlation')
            colormap(gca,posnegcolor)
            ColorbarWithAxis([-0.2 0.2],'rho')
        subplot(2,2,4)
            imagesc(log2(freqs),1:spikes.numcells,spikephasemag(spikephasesort,:))
            xlabel('f (Hz)');ylabel(['Cell - Sorted by ',sortname])
            title('Spike - Phase Coupling')
            LogScale('x',2)
            caxis([0 0.2])
            colorbar

        if figfolder
           NiceSave('SpikeLFPCoupling',figfolder,[]) 
        end

    end
    
    else   %MULTI-CHANNEL FIGURE
        numchans = length(LFP.channels);
        [ratepowerdist,powerbins] = hist(squeeze(ratepowercorr),40);
        
        phasebins = linspace(-pi,pi,40);phasebins = phasebins(1:end-1)+0.5.*diff(phasebins([1 2]));
        [phasehisthist] = hist(squeeze(spikephaseangle),phasebins);
        [phasemaghist,phasemagbins] = hist(squeeze(spikephasemag),40);
        %%
        
    switch nfreqs  
        case 1
    %% Figure: 1 Freq Band
    figure
        subplot(3,2,1)
            %hold on
            imagesc(powerbins,1:numchans,ratepowerdist')
            %axis xy
            hold on
            for pp = 1:numpop
                plot(popcoupling.(pops{pp}).powercorr,1:numchans)
            end
            axis tight
            xlabel('Rate-Power Correlation')
            title('Rate-Power Coupling')
            ylabel('# Cells')
        subplot(3,2,2)
            imagesc(phasebins,1:numchans,phasehisthist')
            hold on
            imagesc(phasebins+2*pi,1:numchans,phasehisthist')
            xlim([-pi 3*pi])
            for pp = 1:numpop
                plot(popcoupling.(pops{pp}).phaseangle,1:numchans,'.')
                plot(popcoupling.(pops{pp}).phaseangle+2*pi,1:numchans,'.')

            end
            title('Spike-Phase Coupling')
        subplot(3,2,4)
            %hold on
            imagesc(phasemagbins,1:numchans,phasemaghist')
            %axis xy
            hold on
            for pp = 1:numpop
                plot(popcoupling.(pops{pp}).powercorr,1:numchans)
            end
            axis tight
            xlabel('Rate-Power Correlation')
            title('Rate-Power Coupling')
            ylabel('# Cells')
            
%         subplot(3,2,3)
%             for pp = 1:numpop
%                 hold on
%                 plot(ratepowercorr(popcellind{pp}),(spikephasemag(popcellind{pp})),'o')
%                 plot(popcoupling.(pops{pp}).powercorr,popcoupling.(pops{pp}).phasemag,'*')
%             end
%             %LogScale('y',10)
%             xlabel('Rate-Power Correlation')
%             ylabel('Spike-Phase Coupling Magnitude')
        
        if figfolder
           NiceSave('SpikeLFPCoupling',figfolder,[]) 
        end

        otherwise
    %% Figure: Spectrum
    posnegcolor = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);

    figure
        subplot(3,2,1)
            hold on
            for pp = 1:numpop
                plot(log2(freqs),popcoupling.(pops{pp}).powercorr)
            end
            plot(log2(freqs([1 end])),[0 0],'k--')
            axis tight
            box off
            LogScale('x',2)
            title('Pop. Synchrony - Power Correlation')
            xlabel('f (Hz)');ylabel('rho')

    %     subplot(3,2,3)
    %         hold on
    %             [hAx,hLine1,hLine2] = plotyy(log2(freqs),cat(1,synchcoupling.phasemag),...
    %                 log2(freqs),mod(cat(1,synchcoupling.phaseangle),2*pi))
    %         LogScale('x',2)
    %         hLine1.LineStyle = 'none';
    %         hLine2.LineStyle = 'none';
    %         hLine1.Marker = 'o';
    %         hLine2.Marker = '.';
    %         title('Pop. Synchrony - Phase Coupling')
    %         xlabel('f (Hz)');
    %         ylabel(hAx(1),'Phase Coupling Magnitude')
    %         ylabel(hAx(2),'Phase Coupling Angle')

        subplot(3,2,3)
            hold on
            for pp = 1:numpop
                plot(log2(freqs),popcoupling.(pops{pp}).phasemag)
            end
            axis tight
            box off
            LogScale('x',2)
            title('Pop. Synchrony - Phase Coupling')
            xlabel('f (Hz)');
            ylabel('Phase Coupling Magnitude')
        subplot(3,2,5)
            hold on
            for pp = 1:numpop
                plot(log2(freqs),popcoupling.(pops{pp}).phaseangle,'o')
                plot(log2(freqs),popcoupling.(pops{pp}).phaseangle+2*pi,'o')
            end
            LogScale('x',2)
            title('Pop. Synchrony - Phase Coupling')
            xlabel('f (Hz)');
            axis tight
            box off
            ylim([-pi 3*pi])
            ylabel('Phase Coupling Angle')


        subplot(2,2,2)
            imagesc(log2(freqs),1:spikes.numcells,ratepowercorr(spikepowersort,:))
            xlabel('f (Hz)');ylabel(['Cell - Sorted by ',sortname])
            LogScale('x',2)
            title('Rate - Power Correlation')
            colormap(gca,posnegcolor)
            ColorbarWithAxis([-0.2 0.2],'rho')
        subplot(2,2,4)
            imagesc(log2(freqs),1:spikes.numcells,spikephasemag(spikephasesort,:))
            xlabel('f (Hz)');ylabel(['Cell - Sorted by ',sortname])
            title('Spike - Phase Coupling')
            LogScale('x',2)
            caxis([0 0.2])
            colorbar

        if figfolder
           NiceSave('SpikeLFPCoupling',figfolder,[]) 
        end

    end
        
    end
end









    %% Spike-Phase Coupling function
    %takes spike times from a single cell and caluclates phase coupling magnitude/angle
    function [phmag,phangle] = spkphase(spkLFP_fn)
        %Spike Times have to be column vector
            if isrow(spkLFP_fn); spkLFP_fn=spkLFP_fn'; end
            if isempty(spkLFP_fn); phmag=nan;phangle=nan; return; end
        %Calculate (power normalized) resultant vector
        rvect = nanmean(abs(spkLFP_fn).*exp(1i.*angle(spkLFP_fn)),1);
        phmag = abs(rvect);
        phangle = angle(rvect);

        %% Example Figure : Phase-Coupling
        % if phmag >0.1
        % figure
        %         rose(phase4spike)
        %        % set(gca,'ytick',[])
        %         hold on
        %          polar(phase4spike,power4spikephase.*40,'k.')
        %         % hold on
        %          %compass([0 phangle],[0 phmag],'r')
        %          compass(rvect.*700,'r')
        %     delete(findall(gcf,'type','text'));
        %     % delete the text objects
        % end        
    end

    %% Power-ISI coupling function
    function [I2] = ISIpowerMutInf(prevISI,nextISI,spkLFP)
        %powerbins = linspace(-1,1.5,8);
        %logISIbins = linspace(-2.5,1,75);
        allISIs = [prevISI;nextISI];
        allLFP = [spkLFP;spkLFP];
        allLFP(isnan(allISIs),:)=[];
        allISIs(isnan(allISIs))=[];
        
        for ff = 1:size(spkLFP,2)

            %ff
%             joint = hist3([log10(allISIs) log2(abs(allLFP(:,ff)))],{logISIbins,powerbins});
%             joint = joint./sum(joint(:));
% 
%             margX = hist(log10(allISIs),logISIbins);
%             margX = margX./sum(margX);
%             margPower = hist(log2(abs(allLFP(:,ff))),powerbins);
%             margPower = margPower./sum(margPower);
%             jointindependent =  bsxfun(@times, margX.', margPower);
% 
% 
%             mutXPow = joint.*log2(joint./jointindependent); % mutual info at each bin
%             totmutXPow(ff) = nansum(mutXPow(:)); % sum of all mutual information 
%             
%            [ I(ff) ] = kernelmi( log10(allISIs)', log2(abs(allLFP(:,ff)))' );
            I2(ff) = mutualinfo(log10(allISIs),log2(abs(allLFP(:,ff))));
            %I3(ff) = mutualinfo((allISIs),(abs(allLFP(:,ff))));
        end
        
        %% Example Figure: ISI-power modulation
%         figure
% %         subplot(2,2,1)
% %             imagesc((logISIbins),powerbins,joint')
% %             axis xy
% %             LogScale('x',10);LogScale('y',2)
% %             xlabel('ISI (s)');ylabel('Power')
% %         subplot(2,2,2)
% %             imagesc((logISIbins),powerbins,jointindependent')
% %             axis xy
% %             %LogScale('xy',10)
% % 
% %             xlabel('ISI (s)');ylabel('Power')
% 
%         subplot(2,2,3)
%             %plot(totmutXPow)
%             hold on
%             %plot(I)
%              plot(I2)
%              %plot(I3)
        
    end
end
