function bz_meta_extracted = bz_database_extract_meta(bz_datasets,j)
% Submit dataset to the Buzsakilab metadata-database
% v1.1
%
% INPUT
% bz_datasets
% structure with subfields named after the mySQL columns. The subfields has
% to have the right format for submitting to the database.
%
%
% By Peter Petersen
% petersen.peter@gmail.com

bz_meta_extracted = bz_datasets;
if ~exist('j')
    j = size(bz_meta_extracted,2);
end
Session = bz_datasets(j).SessionPath;
fprintf(['\n ', num2str(j), '. Animal ' bz_datasets(j).Animal, ', session ', bz_datasets(j).Session, ': '])

% LFP
if exist([Session,'.eeg'])
    bz_meta_extracted(j).LFP = 'Yes';
elseif exist([Session,'.lfp'])
    bz_meta_extracted(j).LFP = 'Yes';
else
    bz_meta_extracted(j).LFP = 'No';
end

% Raw data
if exist(([Session,'.dat']))
    bz_meta_extracted(j).RawData = 'Raw data';
elseif exist(([Session,'.fil']))
    bz_meta_extracted(j).RawData = 'High-pass filtered';
else
    bz_meta_extracted(j).RawData = 'No';
end
if ~strcmp(bz_meta_extracted(j).RawData,'No')
    fprintf(['Data format: ' char(bz_meta_extracted(j).RawData) '. '])
end

% Sampling rate and channel count
if exist(([Session,'.xml']))
    xml = LoadXml(([Session,'.xml']));
    bz_meta_extracted(j).SampleRate = xml.SampleRate;
    bz_meta_extracted(j).lfpSampleRate = xml.lfpSampleRate;
    bz_meta_extracted(j).nChannels = xml.nChannels;
    bz_meta_extracted(j).nElecGps = xml.nElecGps;
    bz_meta_extracted(j).Unitsformat = 'NA';
end

% Number of Cells
bz_meta_extracted(j).Cells = 0;
for i = 1:bz_meta_extracted(j).nElecGps
    if exist(([Session,'.clu.' num2str(i)]))
        cluster_index = load(([Session,'.clu.' num2str(i)]));
        bz_meta_extracted(j).Cells = bz_meta_extracted(j).Cells + sum(unique(cluster_index(2:end))>1);
        bz_meta_extracted(j).Unitsformat = 'Neurosuite';
    end
end
if bz_meta_extracted(j).Cells>0
    fprintf([num2str(bz_meta_extracted(j).Cells) ' units detected. '])
else
    fprintf('No units detected. ')
end

% Duration
if exist(([Session,'.eeg']))
    temp_ = dir(([Session,'.eeg']));
    bz_meta_extracted(j).Duration = temp_.bytes/bz_meta_extracted(j).nChannels/bz_meta_extracted(j).lfpSampleRate/2;
elseif exist(([Session,'.lfp']))
    temp_ = dir(([Session,'.lfp']));
    bz_meta_extracted(j).Duration = temp_.bytes/bz_meta_extracted(j).nChannels/bz_meta_extracted(j).lfpSampleRate/2;
else
    bz_meta_extracted(j).Duration = 0;
end

% Sleep scored
if exist(([Session,'.states.REM'])) | exist(([Session,'.states.SWS'])) | exist(([Session,'.states.Wake']))
    bz_meta_extracted(j).SleepScored = 'Yes';
elseif exist(([Session,'-states.mat']))
    bz_meta_extracted(j).SleepScored = 'Yes';
elseif exist(([Session,'.sts.REM'])) | exist(([Session,'.sts.SWS'])) | exist(([Session,'.sts.Wake']))
    bz_meta_extracted(j).SleepScored = 'Yes';
elseif exist(([Session,'_SleepScore.mat']))
    bz_meta_extracted(j).SleepScored = 'Yes';
elseif exist(([Session,'.SleepState.states.mat']))
    bz_meta_extracted(j).SleepScored = 'Yes';
else
    bz_meta_extracted(j).SleepScored = 'No';
end
if strcmp(bz_meta_extracted(j).SleepScored, 'Yes')
    fprintf('SleepScoring detected. ')
else
    fprintf('No sleepScoring. ')
end

% Ripples
if exist(([Session,'.ripplesALL.event.mat']))
    fprintf('Ripple analysis detected. ')
    bz_meta_extracted(j).RipplesDetected = 'Yes';
else
    fprintf('No ripple analysis detected. ')
    bz_meta_extracted(j).RipplesDetected = 'No';
end

% Tracking
if exist(([Session,'.whl']))
    fprintf('whl tracking files detected. ')
    bz_meta_extracted(j).Tracking = 'Camera';
    bz_meta_extracted(j).TrackingFormat = 'whl';
    bz_meta_extracted(j).TrackingFrameRate = 39.06;
else
    bz_meta_extracted(j).Tracking = 'No';
    bz_meta_extracted(j).TrackingFormat = nan;
end

% Buzcode
if exist(([Session,'.sessionInfo.mat']))
    fprintf('Buzcode sessionInfo file detected. ')
    sessionInfo = load([Session,'.sessionInfo.mat']);
    if isfield(sessionInfo.sessionInfo,'region') && ~isempty(sessionInfo.sessionInfo.region)
        bz_meta_extracted(j).BrainRegion = strjoin(unique(sessionInfo.sessionInfo.region),',');
    end
    if isfield(sessionInfo.sessionInfo,'depth') && ~isempty(sessionInfo.sessionInfo.region)
        bz_meta_extracted(j).Depth = char(sessionInfo.sessionInfo.depth);
    end
    if isfield(sessionInfo.sessionInfo,'Date')
        bz_meta_extracted(j).Date = sessionInfo.sessionInfo.Date;
    end
    bz_meta_extracted(j).FileFormat = 'Buzcode';
    clear sessionInfo
end

if exist(([Session,'.behavior.mat']))
    fprintf('Buzcode behavior file detected. ')
    behavior = load([Session,'.behavior.mat']);
    bz_meta_extracted(j).TrackingFormat = 'buzcode';
    if isfield(behavior.behavior,'acquisitionsystem') && ~isempty(behavior.behavior.acquisitionsystem)
        bz_meta_extracted(j).Tracking = behavior.behavior.acquisitionsystem;
    end
    if isfield(behavior.behavior,'samplingRate') && ~isempty(behavior.behavior.samplingRate)
        bz_meta_extracted(j).TrackingFrameRate = behavior.behavior.samplingRate;
    end
    if isfield(behavior.behavior,'events') && isfield(behavior.behavior.events,'trials') && ~isempty(behavior.behavior.events.trials)
    bz_meta_extracted(j).BehavioralTrials = size(behavior.behavior.events.trials,2);
    end
    clear behavior
end

% Brendon's standard
if exist(([Session,'_BasicMetaData.mat']))
    fprintf('_BasicMetaData file detected (Brendons format). ')
    BasicMetaData = load([Session,'_BasicMetaData.mat']);
    if isfield(BasicMetaData,'CortexRegion')
        bz_meta_extracted(j).BrainRegion = strjoin(unique(BasicMetaData.CortexRegion),',');
    end
    if isfield(BasicMetaData,'bmd') && isfield(BasicMetaData.bmd,'Par') && isfield(BasicMetaData.bmd.Par,'Date')
        bz_meta_extracted(j).Date = BasicMetaData.bmd.Par.Date;
    end
    bz_meta_extracted(j).FileFormat = 'Brendon_metadata';
    clear bmd CortexRegion
end

% Neurodata Without Borders
