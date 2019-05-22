function sessionInfo = get_SessionInfoFromNWB(nwb2)

%% get_SessionInfoFromNWB populates a variable similar to sessionInfo standards.
% Not all fields are filled (most are commented out). 
% Only those that were needed so far in the functions tested. 
% More can be filled on demand.




%% Get the type of the recording (its key will be used to get info from the nwb file)
% There will be 2 logical values:
% RawDataPresent
% LFPDataPresent

% The keys are not finalized yet from the NWB format so some updating will
% be needed.


%% Populate fields


% First check if there are raw data present

try
    all_raw_keys = keys(nwb2.acquisition);

    if ~isempty(all_raw_keys)
        for iKey = 1:length(all_raw_keys)
            if ismember(all_raw_keys{iKey}, {'ECoG','bla bla bla'})   %%%%%%%% ADD MORE HERE, DON'T KNOW WHAT THE STANDARD FORMATS ARE
                iRawDataKey = iKey;
                RawDataPresent = 1;
            else
                RawDataPresent = 0;
            end
        end
        % nwb2.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get('all_lfp').data
    else
        RawDataPresent = 0;
    end
catch
    RawDataPresent = 0;
end




try
    % Check if the data is in LFP format
    all_lfp_keys = keys(nwb2.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries);

    if ~isempty(all_lfp_keys)
        for iKey = 1:length(all_lfp_keys)
            if ismember(all_lfp_keys{iKey}, {'lfp','bla bla bla'})   %%%%%%%% ADD MORE HERE, DON'T KNOW WHAT THE STANDARD FORMATS ARE
                iLFPDataKey = iKey;
                LFPDataPresent = 1;
                break % Once you find the data don't look for other keys/trouble
            else
                LFPDataPresent = 0;
            end
        end
    else
        LFPDataPresent = 0;
    end
catch
    LFPDataPresent = 0;
end


if ~RawDataPresent && ~LFPDataPresent
    error 'There is no data in this .nwb - Maybe check if the Keys are labeled correctly'
end



%% Creare SessionInfo
sessionInfo = struct;

if RawDataPresent
    sessionInfo.nChannels      = nwb2.acquisition.get(all_raw_keys{iRawDataKey}).data.dims(2);
    sessionInfo.samples_NWB    = nwb2.acquisition.get(all_raw_keys{iRawDataKey}).data.dims(1);
    sessionInfo.rates.wideband = nwb2.acquisition.get(all_raw_keys{iRawDataKey}).starting_time_rate; % THIS SHOULD BE 20,000 HZ FOR THE TUTORIAL
    sessionInfo.rates.lfp      = 1250;    %1250 -  DEFAULT - CHECK THIS
    sessionInfo.lfpSampleRate  = 1250;    %1250 -  DEFAULT - CHECK THIS % tH lfpSampleRate bypasses 

elseif LFPDataPresent
    sessionInfo.nChannels      = nwb2.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get(all_lfp_keys{iLFPDataKey}).data.dims(2);
    sessionInfo.samples_NWB    = nwb2.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get(all_lfp_keys{iLFPDataKey}).data.dims(1);
    sessionInfo.rates.wideband = nwb2.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get(all_lfp_keys{iLFPDataKey}).starting_time_rate;
    sessionInfo.rates.lfp      = nwb2.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get(all_lfp_keys{iLFPDataKey}).starting_time_rate;    %I assign the LFP sampling rate that was already used. Not sure yet if a different value than 1250 Hz will cause problems
    sessionInfo.lfpSampleRate  = nwb2.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get(all_lfp_keys{iLFPDataKey}).starting_time_rate;    %I assign the LFP sampling rate that was already used. Not sure yet if a different value than 1250 Hz will cause problems
    
    sessionInfo.LFPDataPresent = LFPDataPresent;            % This is added for easy retrieval from bz_GetLFP
    sessionInfo.LFPDataKey     = all_lfp_keys{iLFPDataKey}; % This is added for easy retrieval from bz_GetLFP
end



% Get the ElectrodeGroups
% There is a lot of redundancy on this one

% electrode_description = nwb2.general_extracellular_ephys_electrodes.vectordata.get('electrode_description').data;
% filtering = nwb2.general_extracellular_ephys_electrodes.vectordata.get('filtering').data;
% group = nwb2.general_extracellular_ephys_electrodes.vectordata.get('group').data;
% group_name = nwb2.general_extracellular_ephys_electrodes.vectordata.get('group_name').data;
% imp = nwb2.general_extracellular_ephys_electrodes.vectordata.get('imp').data.load;
location = nwb2.general_extracellular_ephys_electrodes.vectordata.get('location').data.load;
% shank = nwb2.general_extracellular_ephys_electrodes.vectordata.get('shank').data.load;
% x = nwb2.general_extracellular_ephys_electrodes.vectordata.get('x').data.load;
% y = nwb2.general_extracellular_ephys_electrodes.vectordata.get('y').data.load;
% z = nwb2.general_extracellular_ephys_electrodes.vectordata.get('z').data.load;



% nGroups = sum(unique(shank)>0); % -1 doesn't belong to a Shank Group
nGroups          = length(unique(nwb2.general_extracellular_ephys_electrodes.vectordata.get('group_name').data.load));
uniqueGroupNames = unique(nwb2.general_extracellular_ephys_electrodes.vectordata.get('group_name').data.load);
groupNames       = nwb2.general_extracellular_ephys_electrodes.vectordata.get('group_name').data.load;


sessionInfo.spikeGroups.nGroups  = nGroups;
sessionInfo.spikeGroups.nSamples = ones(1,sessionInfo.spikeGroups.nGroups)*32; % The file I found had 32 here


id = nwb2.general_extracellular_ephys_electrodes.vectordata.get('amp_channel').data.load;


sessionInfo.spikeGroups.groups = cell(1,sessionInfo.spikeGroups.nGroups);
for iGroup = 1:sessionInfo.spikeGroups.nGroups
    sessionInfo.spikeGroups.groups{iGroup} = id(ismember(groupNames, uniqueGroupNames{iGroup}))';
    
    % Redundant
    for iChannel = 1:length(id(ismember(groupNames, uniqueGroupNames{iGroup}))')
        sessionInfo.ElecGp{1,iGroup}.channel{1,iChannel} = sessionInfo.spikeGroups.groups{iGroup}(iChannel);
    end
    
    % Redundant
    sessionInfo.SpkGrps(iGroup).Channels   = id(ismember(groupNames, uniqueGroupNames{iGroup}))';
    sessionInfo.SpkGrps(iGroup).nSamples   = 32;
    sessionInfo.SpkGrps(iGroup).PeakSample = 16; 
    sessionInfo.SpkGrps(iGroup).nFeatures  = 3; 
    
    
end


% Add Region Info
% If the recording has additional channels (behavioral channels - they wouldn't have a location field - 
% I assume here that the numbering of the channels starts with the electrophysiological, and the behavioral are concatenated after)
% Thats why I use length(groupNames) here and not sessionInfo.nChannels

for iChannel = 1:length(groupNames) 
    sessionInfo.region{iChannel} = location{iChannel};
end


% Get the rest of the Info
sessionInfo.nBits          = 16; % ASSUMING THAT NWB HAS SAVED DATA IN INT16 PRECISION
sessionInfo.rates.video    = 0;
sessionInfo.FileName       = nwb2.identifier;%%%%%%%%%%%%%%%% no extension - I DON'T USE THE FILENAME OF THE NWB HERE JUST IN CASE SOMEONE CHANGED IT. 
% sessionInfo.SampleTime     = 50; % 50 no idea
sessionInfo.nElecGps       = nGroups; % 13
% sessionInfo.ElecGp         = []; % 1x13 cell (struct with 1x12 cell inside)
% sessionInfo.HiPassFreq     = ;% probably the one from the LFP conversion
% sessionInfo.Date           = 
% sessionInfo.VoltageRange   = % 20
% sessionInfo.Amplification  = % 1000
% sessionInfo.Offset         = 0;
% sessionInfo.AnatGrps       =
% sessionInfo.spikeGroups.groups        = {1:32};
% sessionInfo.spikeGroups.nSamples = 1; % I ADDED THIS FOR bz_GetSpikes
sessionInfo.channels       =  0:sessionInfo.nChannels-1 ;% 1x128 % starts from 0              THIS IS USED IN THE bz_GetLFP in the end to assign channels to regions
% sessionInfo.lfpChans       =
% sessionInfo.thetaChans     =
% sessionInfo.region         = % cell 1x128
% sessionInfo.depth          = 1;   % 3324 - single value
% sessionInfo.ca1            = 116; % 116 - single value
% sessionInfo.ca3            = [] % []
% sessionInfo.ls             = [];
% sessionInfo.animal         = 'MONKEY'% string
% sessionInfo.refChan        = 112; % single value

sessionInfo.useRaw           = RawDataPresent; % I add this to choose which data to convert to .lfp files



end

