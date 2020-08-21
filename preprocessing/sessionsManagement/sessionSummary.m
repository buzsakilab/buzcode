
% sessionSummary.m

% This scritp runs preliminary descriptive analysis to get a summary overview
% of the session.

% HISTORY:
% - Based on Manu Valero-BuzsakiLab 2019
% - Some reorganization, pre-processing steps transferred to expPipeline: 5/20, AntonioFR


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mkdir('SummaryFigures'); % create folder
close all

% manually adjust analog pulse detection
if ~isempty(analogCh)
	[pulses] = bz_getAnalogPulses('analogCh',analogCh,'manualThr',false);
else
    pulses.intsPeriods = [];
end

% get SleepStateScoring
targetFile = dir('*SleepState.states.mat');
if isempty(targetFile)
    try 
        SleepScoreMaster(pwd,'noPrompts',true,'ignoretime',pulses.intsPeriods); % try to sleep score
    catch
        disp('Problem with SleepScore skyping...');
    end
end

%% 1. Spike-waveform, ACG, spatial position and CCG
try
    disp('Spike-waveform, ACG and cluster location...');
    % spikes = loadSpikes;
    spikes = bz_LoadPhy('noPrompts',true,'verbose',true,'nWaveforms',200);

    % plot spikes summary
    disp('Plotting spikes summary...');
    if exist('chanMap.mat','file')
        load('chanMap.mat','xcoords','ycoords');
        xcoords = xcoords - min(xcoords); xcoords = xcoords/max(xcoords);
        ycoords = ycoords - min(ycoords); ycoords = ycoords/max(ycoords); 
    else
        xcoords = NaN;
        ycoords = NaN;
    end
    figure
    for jj = 1:size(spikes.UID,2)
        fprintf(' **Spikes from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
        dur = 0.08;
        set(gcf,'Position',[100 -100 2500 1200])
        subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
        ccg=CCG(spikes.times,[],'binSize',0.001,'duration',0.08);
        xt = linspace(-dur/2*1000,dur/2*1000,size(ccg,1));
        area(xt,ccg(:,jj,jj),'LineStyle','none');
        try 
            xt2 = linspace(dur/3*1000+5,dur/2*1000+5+80,size(spikes.filtWaveform{jj},1)); % waveform
            spk = spikes.filtWaveform{jj} - min(spikes.filtWaveform{jj}); spk = spk/max(spk) * max(ccg(:,jj,jj));
            hold on
            plot(xt2,spk)

            plot((xcoords*30) + dur/2*1000+5+60, ycoords*max(ccg(:,jj,jj)),'.','color',[.8 .8 .8],'MarkerSize',5); % plotting xml
            plot((xcoords(find(spikes.maxWaveformCh(jj)==sessionInfo.channels))*30) + dur/2*1000+5+60,...
                ycoords(find(spikes.maxWaveformCh(jj)==sessionInfo.channels))*max(ccg(:,jj,jj)),'.','color',[.1 .1 .1],'MarkerSize',10); % plotting xml
        end
        title(num2str(jj),'FontWeight','normal','FontSize',10);
        if jj == 1
            ylabel('Counts/ norm');
        elseif jj == size(spikes.UID,2)
            xlabel('Time (100 ms /1.5ms)');
        else
            set(gca,'YTick',[],'XTick',[]);
        end
    end
    saveas(gcf,'SummaryFigures\spikes.png');

    win = [-0.3 0.3];
    disp('Plotting CCG...');
    figure;
    set(gcf,'Position',[100 -100 2500 1200])
    [allCcg, t] = CCG(spikes.times,[],'binSize',0.005,'duration',0.6);
    indCell = [1:size(allCcg,2)];
    for jj = 1:size(spikes.UID,2)
        fprintf(' **CCG from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
        subplot(7,ceil(size(spikes.UID,2)/7),jj);
        cc = zscore(squeeze(allCcg(:,jj,indCell(indCell~=jj)))',[],2); % get crosscorr
        imagesc(t,1:max(indCell)-1,cc)
        set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
        hold on
        zmean = mean(zscore(squeeze(allCcg(:,jj,indCell(indCell~=jj)))',[],2));
        zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
        plot(t, zmean,'k','LineWidth',2);
        xlim([win(1) win(2)]); ylim([0 max(indCell)-1]);
        title(num2str(jj),'FontWeight','normal','FontSize',10);

        if jj == 1
            ylabel('Cell');
        elseif jj == size(spikes.UID,2)
            xlabel('Time (s)');
        else
            set(gca,'YTick',[],'XTick',[]);
        end
    end
    saveas(gcf,'SummaryFigures\CrossCorr.png'); 
catch
    warning('Error on Spike-waveform, autocorrelogram and cluster location! ');
end


%% 2. Psth and CSD from analog-in inputs
try 
    disp('Psth and CSD from analog-in inputs...');
    % copy NPY file to Kilosort folder, for phy psth option
    writeNPY(pulses.timestamps(:,1),'pulTime.npy'); % if error, transpose!
    kilosort_path = dir('*Kilosort*');
    try copyfile('pulTime.npy', strcat(kilosort_path.name,'\','pulTime.npy')); end% copy pulTime to kilosort folder
    xml = LoadParameters;

    eventIDList = unique(pulses.eventID);
    if exist('pulses') && ~isempty(pulses.timestamps(:,1))
        for mm = 1:length(eventIDList)
            fprintf('Stimulus %3.i of %3.i \n',mm, length(eventIDList)); %\n
            st = pulses.timestamps(pulses.eventID == eventIDList(mm),1);
            % CSD
            figure
            for jj = 1:size(xml.AnatGrps,2)
                lfp = bz_GetLFP(xml.AnatGrps(jj).Channels(1:numel(xml.AnatGrps(jj).Channels)-mod(numel(xml.AnatGrps(jj).Channels),8)),'noPrompts', true);
                twin = 0.02;
                [csd,lfpAvg] = bz_eventCSD(lfp,st,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
                taxis = linspace(-twin,twin,size(csd.data,1));
                cmax = max(max(csd.data)); 
                subplot(1,size(xml.AnatGrps,2),jj);
                contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
                set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('STIMULATION, Shank #',num2str(jj)),'FontWeight','normal'); 
                colormap jet; try caxis([-cmax cmax]); end
                hold on
                for kk = 1:size(lfpAvg.data,2)
                    plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k')
                end
            end
            saveas(gcf,['SummaryFigures\stimulation_',num2str(mm), '.png']);

            % PSTH
            win = [-0.1 0.5];
            disp('Plotting spikes raster and psth...');
            figure;
            set(gcf,'Position',[100 -100 2500 1200])
            for jj = 1:size(spikes.UID,2)
                fprintf(' **Pulses from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
                rast_x = []; rast_y = [];
                for kk = 1:length(st)
                    temp_rast = spikes.times{jj} - st(kk);
                    temp_rast = temp_rast(temp_rast>win(1) & temp_rast<win(2));
                    rast_x = [rast_x temp_rast'];
                    rast_y = [rast_y kk*ones(size(temp_rast))'];
                end
                [stccg, t] = CCG({spikes.times{jj} st},[],'binSize',0.005,'duration',1);
                subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                plot(rast_x, rast_y,'.','MarkerSize',1)
                hold on
                plot(t(t>win(1) & t<win(2)), stccg(t>win(1) & t<win(2),2,1) * kk/max(stccg(:,2,1))/2,'k','LineWidth',2);
                xlim([win(1) win(2)]); ylim([0 kk]);
                title(num2str(jj),'FontWeight','normal','FontSize',10);

                if jj == 1
                    ylabel('Trial');
                elseif jj == size(spikes.UID,2)
                    xlabel('Time (s)');
                else
                    set(gca,'YTick',[],'XTick',[]);
                end
            end
            saveas(gcf,['SummaryFigures\AnalogPulse_',num2str(mm) ,'.png']); 
        end  
    end
catch
    warning('Error on Psth and CSD from analog-in inputs! ');
end


%% 3. Slow-waves CSD and PSTH
try 
    disp('Slow-waves CSD and PSTH...');
    UDStates = detectUD;
    xml = LoadParameters;

    % CSD
    twin = 0.2;
    figure
    for jj = 1:size(xml.AnatGrps,2)
        lfp = bz_GetLFP(xml.AnatGrps(jj).Channels(1:numel(xml.AnatGrps(jj).Channels)-mod(numel(xml.AnatGrps(jj).Channels),8)),'noPrompts', true);
        evs = UDStates.timestamps.DOWN(:,1); evs(evs>lfp.duration);
        [csd,lfpAvg] = bz_eventCSD(lfp,evs,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
        taxis = linspace(-0.2,0.2,size(csd.data,1));
        cmax = max(max(csd.data)); 
        subplot(1,size(xml.AnatGrps,2),jj);
        contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
        set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('DOWN-UP, Shank #',num2str(jj)),'FontWeight','normal'); 
        colormap jet; caxis([-cmax cmax]);
        hold on
        for kk = 1:size(lfpAvg.data,2)
            plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k')
        end
    end
    saveas(gcf,'SummaryFigures\downUp.png');

    % PSTH
    st = UDStates.timestamps.DOWN;
    win = [-0.2 0.2];
    figure
    set(gcf,'Position',[100 -100 2500 1200])
    for jj = 1:size(spikes.UID,2)
        fprintf(' **UD from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
        rast_x = []; rast_y = [];
        for kk = 1:length(st)
            temp_rast = spikes.times{jj} - st(kk);
            temp_rast = temp_rast(temp_rast>win(1) & temp_rast<win(2));
            rast_x = [rast_x temp_rast'];
            rast_y = [rast_y kk*ones(size(temp_rast))'];
        end
        [stccg, t] = CCG({spikes.times{jj} st},[],'binSize',0.005,'duration',1);
        subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
        plot(rast_x, rast_y,'.','MarkerSize',1)
        hold on
        plot(t(t>win(1) & t<win(2)), stccg(t>win(1) & t<win(2),2,1) * kk/max(stccg(:,2,1))/2,'k','LineWidth',2);
        xlim([win(1) win(2)]); ylim([0 kk]);
        title(num2str(jj),'FontWeight','normal','FontSize',10);

        if jj == 1
            ylabel('Trial');
        elseif jj == size(spikes.UID,2)
            xlabel('Time (s)');
        else
            set(gca,'YTick',[],'XTick',[]);
        end
    end
    saveas(gcf,'SummaryFigures\psthUDStates.png'); 
catch
    warning('Error on Psth and CSD from analog-in inputs! ');
end


%% 6. BEHAVIOUR
try 
    getSessionTracking('convFact',0.1149,'roiTracking','manual'); 
    getSessionArmChoice;
    getSessionLinearize;  
catch
    warning('It has not been possible to run the behaviour code...');
end


