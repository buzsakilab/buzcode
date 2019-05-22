%% Example conversion to NWB of the YutaMouse41-150903 dataset

% Folder that contains all the files - This is the only input needed
folder_path = 'F:\NWBtoBuzcode\YutaMouse41-150903';
% folder_path = 'C:\Users\McGill\Documents\GitHub\buzcode\tutorials\exampleDataStructs\20170505_396um_0um_merge';



%% Start Adding fields to NWB

% Get info from the xml file
xml = GetXMLInfo(folder_path);

% Add general info to the NWB file
nwb = GeneralInfo(xml);

% Add electrode info
nwb = addElectrodeInfo(xml, nwb);

% Add units info - By default, the spike waveforms are added to the file
nwb = addUnitsInfo_Neuroscope(xml,nwb);

% Add stimulation events
nwb = addEvents_Yuta(xml,nwb);

% Add behavioral info/channels
nwb = addBehavior_Yuta(xml,nwb);

% Add electrophysiological channels
nwb = addElectrophysiology(xml, nwb);

% Add epochs
nwb = addEpochs_Yuta(xml,nwb);

% Add trials
nwb = addTrials_Yuta(xml,nwb);

% Add channels based on the Yuta spreadsheet
nwb = addSpecial_YutaMouse_recordings(xml,nwb);


%% Export to nwb
% nwbExport(nwb, 'YutaMouse41_converted.nwb')
nwbExport(nwb, 'F:\NWBtoBuzcode\YutaMouse41-150903\YutaMouse41.nwb')

