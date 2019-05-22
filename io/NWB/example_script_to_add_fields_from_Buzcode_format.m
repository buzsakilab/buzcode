

% test standard formats
folder_path = 'C:\Users\McGill\Documents\GitHub\buzcode\tutorials\exampleDataStructs\20170505_396um_0um_merge';


%% Start Adding fields to NWB

% Get info from the rhd file
xml = GetXMLInfo(folder_path);

% Add general info to the NWB file
nwb = GeneralInfo(xml);

% Add electrode information
nwb = addElectrodeInfo(xml, nwb);

% Add Behavior information
nwb = addBehavior(xml, nwb);

% Add Behavior information
nwb = addTrials(xml, nwb);

% Add events
nwb = addEvents(xml, nwb);

% Add electrophysiological data (.lfp, .eeg)
nwb = addElectrophysiology(xml, nwb);

% Add units info
nwb = addUnitsInfo(xml, nwb);

%% Export to nwb
nwbExport(nwb, 'F:\NWBtoBuzcode\test_Buzcode_Standards.nwb')








