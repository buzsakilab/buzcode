%% Example conversion to NWB of the YutaMouse41-150903 dataset

% Folder that contains all the files - This is the only input needed
folder_path = 'F:\NWBtoBuzcode\YutaMouse41-150903';

%% Start Adding fields to NWB

% Get info from the xml file
xml = Neuroscope2NWB.GetXMLInfo(folder_path);

% Add general info to the NWB file
nwb = Neuroscope2NWB.GeneralInfo(xml);

% Add electrode info
nwb = Neuroscope2NWB.getElectrodeInfo(xml, nwb);

% Add units info - By default, the spike waveforms are added to the file
nwb = Neuroscope2NWB.getUnitsInfo(xml,nwb);

% Add stimulation events
nwb = Neuroscope2NWB.getEvents(xml,nwb);

% Add behavioral info/channels
nwb = Neuroscope2NWB.getBehavior(xml,nwb);

% Add electrophysiological channels
nwb = Neuroscope2NWB.getElectrophysiology(xml, nwb);

% Add epochs
nwb = Neuroscope2NWB.getEpochs(xml,nwb);

% Add trials
nwb = Neuroscope2NWB.getTrials(xml,nwb);

% Add channels based on the Yuta spreadsheet
nwb = Neuroscope2NWB.special_YutaMouse_recordings(xml,nwb);


%% Export to nwb
nwbExport(nwb, 'YutaMouse41_converted.nwb')

