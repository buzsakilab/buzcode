% this script is a brief tutorial for mapping the firing rates of neurons
% onto the behavioral positions occupied by an animal
%
% This script assumes that you are currently in the
% /buzcode/tutorials/exampleDataStructs/20170505_396um_0um_merge folder, 
% it will not work otherwise!

path = strsplit(pwd,filesep);
if ~strcmp(path{end},'tutorials') & ~strcmp(path{end-1},'buzcode') 
   error('this script assumes you are running it from the /buzcode/tutorials folder...') 
end
cd(['exampleDataStructs' filesep '20170505_396um_0um_merge'])

%% loading data
sessionInfo = bz_getSessionInfo;
lfpChan = sessionInfo.ca1;
spikes = bz_GetSpikes('noprompts',true);
try
    lfp = bz_GetLFP(lfpChan);
catch
    error('couldnt load LFP, have you run the "download_DATA" script yet?')
end

behavior = bz_LoadBehavior(pwd,'track');
% what does the behavior look like?
subplot(3,2,1)
bz_plotTrials(behavior,'condition',1)

%% mapping spiking onto behavior
[firingMaps] = bz_firingMap1D(spikes,behavior,4,'savemat',false);
% let's look at some ratemaps...
subplot(3,2,2)
imagesc(squeeze(firingMaps.rateMaps{1}(96,:,:)))
ylabel('trial #')
xlabel('position')

[phaseMaps] = bz_phaseMap1D(spikes,behavior,lfp,4,'savemat',false);
% phase precession, eh?
subplot(3,2,4)
scatter(phaseMaps.phaseMaps{1}{96}(:,1),phaseMaps.phaseMaps{1}{96}(:,end),'.k')
hold on
scatter(phaseMaps.phaseMaps{1}{96}(:,1),phaseMaps.phaseMaps{1}{96}(:,end)+2*pi,'.k')
axis([0 200 -pi pi*3])
ylabel('theta phase')
xlabel('position')


subplot(3,2,3)
% pretty phasemap for rachel.... coming when someone has the motivation :)

%% detecting place fields
% should we find some place fields?
for condition = 1:length(firingMaps.rateMaps)
     fields{condition} = bz_getPlaceFields1D(firingMaps.rateMaps{condition},'minPeakRate',2,'percentThreshold',.15);
end

% typing 'fields{1}{96}{1}' will give you an idea of the type of data
% stored for each place field
fieldStart = fields{1}{96}{1}.start;
fieldStop = fields{1}{96}{1}.stop;
fieldCenterMass = fields{1}{96}{1}.COM;
% fieldStart = fields{1}{96}{1}.start;


subplot(3,2,3)
plot(squeeze(mean(firingMaps.rateMaps{1}(96,:,:),2)),'k')
line([fieldStart fieldStart],[0 25],'color','g')
line([fieldStop fieldStop],[0 25],'color','g')
line([fieldCenterMass fieldCenterMass],[0 25],'color','r')
title('G - start/stop, R - Center of mass')
ylabel('avg. firing rate')
xlabel('position')

