

%% Test BUZSAKI FUNCTIONS and NWB FUNCTIONS

% For the NWB functions, the nwb file_path should be added as an input
% (no need to have only one .nwb in the same folder)
nwb_file = 'C:\Users\McGill\Documents\GitHub\matnwb\Nas\YutaMouse41\YutaMouse41\YutaMouse41.nwb';


% If the nwb_file is not specified, the functions will search for a unique
% .nwb within the current directory

% The current Buzcode functions and the NWB functions are differentiated by
% the _NWB suffix.


% Konstantinos Nasiotis 2019


%% GET_LFP

% channel ID 5 (= # 6)
tic
lfp = bz_GetLFP(5);
disp(['Loading a single channel with Buzsaki code takes: ' num2str(toc)])   % 18.0 seconds for a 40369 seconds recording at 1250 Hz

tic
lfp_NWB = bz_GetLFP_NWB(5, 'nwb_file', nwb_file);
disp(['Loading a single channel with NWB code takes: ' num2str(toc)])       % 4.3 seconds for a 40369 seconds recording at 1250 Hz

figure(1);plot(lfp.data)
figure(2);plot(lfp_NWB.data)


% channel ID 5 and 10 (= #6 and #11)
tic
lfp = bz_GetLFP([5 10]);
disp(['Loading a single channel with Buzsaki code takes: ' num2str(toc)])   % 18 seconds for a 40369 seconds recording at 1250 Hz

tic
lfp_NWB = bz_GetLFP_NWB([5 10], 'nwb_file', nwb_file);
disp(['Loading a single channel with NWB code takes: ' num2str(toc)])       % 4.3 seconds for a 40369 seconds recording at 1250 Hz

figure(1);subplot(1,2,1);plot(lfp.data(:,1)); subplot(1,2,2);plot(lfp.data(:,2))
figure(2);subplot(1,2,1);plot(lfp_NWB.data(:,1)); subplot(1,2,2);plot(lfp_NWB.data(:,2))



% channel ID 5 (= # 6), from 0 to 120 seconds
lfp     = bz_GetLFP(5,'restrict',[0 120]);
lfp_NWB = bz_GetLFP_NWB(5,'restrict',[0 120]);

figure(1);plot(lfp.data)
figure(2);plot(lfp_NWB.data)



% same, plus from 240.2 to 265.23 seconds
lfp     = bz_GetLFP(5,'restrict',[0 120;240.2 265.23]);
lfp_NWB = bz_GetLFP_NWB(5, 'restrict',[0 120;240.2 265.23],'nwb_file', nwb_file);


figure(1);subplot(1,2,1); plot(lfp(1).data) ; subplot(1,2,2); plot(lfp(2).data)
figure(2);subplot(1,2,1); plot(lfp_NWB(1).data) ; subplot(1,2,2); plot(lfp_NWB(2).data)



% Downsample recording
lfp     = bz_GetLFP    (5, 'restrict',[0 120], 'downsample', 2);
lfp_NWB = bz_GetLFP_NWB(5, 'restrict',[0 120], 'downsample', 2);

figure(1);plot(lfp.data)
figure(2);plot(lfp_NWB.data)



%% GET SPIKES

% Select ALL Neurons
spikes     = bz_GetSpikes('UID',[2 4 7]);
spikes_NWB = bz_GetSpikes_NWB('UID',[2 4 7]);
spikes_NWB = bz_GetSpikes_NWB('nwb_file', nwb_file, 'UID',[2 4 7]);

% Select Neurons: #2, #4, #7
spikes     = bz_GetSpikes('UID',[2 4 7]);
spikes_NWB = bz_GetSpikes_NWB('nwb_file', nwb_file, 'UID',[2 4 7]);

% Select Neurons from spikeGroups: #2, #4
spikes     = bz_GetSpikes('spikeGroups', [2,4]);
spikes_NWB = bz_GetSpikes_NWB('nwb_file', nwb_file, 'spikeGroups', [2,4]);

% Select Neurons from region: 'unknown'
spikes     = bz_GetSpikes('region', 'unknown');
spikes_NWB = bz_GetSpikes_NWB('nwb_file', nwb_file, 'region', 'unknown');

% Loads and saves all spiking data into buzcode format .cellinfo. struct
spikes = bz_GetSpikes('saveMat',true); 
spikes_NWB = bz_GetSpikes_NWB('nwb_file', nwb_file, 'saveMat',true);



%% Load behavior

% Test that it works
behavior = bz_LoadBehavior;
behavior = bz_LoadBehavior_NWB;
behavior = bz_LoadBehavior_NWB('nwb_file', nwb_file);
behavior = bz_LoadBehavior_NWB('nwb_file', nwb_file, 'EightMazePosition_norm_spatial_series');
behavior = bz_LoadBehavior_NWB('EightMazePosition_norm_spatial_series');


%% Load events
basePath = 'C:\Users\McGill\Documents\GitHub\matnwb\Nas\YutaMouse41\YutaMouse41';
events = bz_LoadEvents(basePath, 'YutaMouse41.PulseStim_0V_10021ms_LD0.events');


events = bz_LoadEvents_NWB;
events = bz_LoadEvents_NWB('nwb_file', nwb_file);
events = bz_LoadEvents_NWB('nwb_file', nwb_file, 'PulseStim_5V_77777ms_LD12');
events = bz_LoadEvents_NWB('PulseStim_5V_77777ms_LD12');


%% bz_firingMap1D

[firingMaps] = bz_firingMap1D(spikes_NWB,behavior,4,'savemat',false);
% let's look at some ratemaps...
subplot(3,2,2)
imagesc(squeeze(firingMaps.rateMaps{1}(96,:,:)))
ylabel('trial #')
xlabel('position')
