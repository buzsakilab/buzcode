function DatsMetadata = bz_DatFileMetadata(basepath)
% Store basic info about the original dat files: names and bytes
% Brendon Watson 2016-8

%% Input and directory handling 
if ~exist('basepath','var')
    basepath = cd;
elseif isempty(basepath)
    basepath = cd;
end

basename = bz_BasenameFromBasepath(basepath);


%% First find out what system was used to record this... based on files present
% check if amplipex (.metas present) or intan (rhd.infos present)... or neither
recsys = '';

%Amplipex case
dmeta = dir(fullfile(basepath,'*.meta'));
if ~isempty(dmeta)
    recsys = 'Amplipex';
end

%Intan case
drhd = dir(basepath);
for didx = length(drhd):-1:1
   ok = 0;
   if drhd(didx).isdir
       t = dir(fullfile(basepath,drhd(didx).name,'info.rhd'));
       if ~isempty(t)%only use folders with info.rhd inside
            recsys = 'Intan';
            ok = 1;
       end
   end
   if ~ok
        drhd(didx) = [];
   end
end

%% Now gather info
switch recsys
    case 'Amplipex'%if from an amplipex
        d = dmeta;
        DatsMetadata.RecordingSystem = 'Amplipex';
        for idx = 1:length(d);
            DatsMetadata.Recordings.Names{idx} = d(idx).name(1:end-5);
            tdat = strcat(d(idx).name(1:end-5),'.dat');
            t = dir(fullfile(basepath,tdat));
            if ~isempty(t) %if original .dat still around
                DatsMetadata.Recordings.Bytes(idx) = t.bytes;
            else%if no .dat around (ie anymore), get bytes from .meta 
                disp([tdat ' not found']);
                DatsMetadata.Recordings.Bytes{idx} = str2num(bz_ReadAmplipexMetafileAspects(fullfile(basepath,d(idx).name),'filebytes'));
            end
            NChannels(idx) = str2num(bz_ReadAmplipexMetafileAspects(fullfile(basepath,d(idx).name),'NumChans'));
            SamplingRate(idx) = str2num(bz_ReadAmplipexMetafileAspects(fullfile(basepath,d(idx).name),'SamplingRate'));
            Amplification(idx) = str2num(bz_ReadAmplipexMetafileAspects(fullfile(basepath,d(idx).name),'Gain'));
        end
        if ~isempty(NChannels) %if some files were found
            %Channels
            if sum(abs(diff(NChannels)))%if different numbers of channels differed across recordings in session, error out
                warning('Rercordings contain different numbers of channels - per .meta files.  Will use count from first recording.')
            end
            DatsMetadata.Parameters.NumberOfChannels = NChannels(1);

            %Sampling Rate
            if sum(abs(diff(SamplingRate)))%if different numbers of channels differed across recordings in session, error out
                warning('Rercordings contain different SamplingRates - per .meta files.  Will use number from first recording.')
            end
            DatsMetadata.Parameters.SampleRate = SamplingRate(1);

            %Amplification
            if sum(abs(diff(Amplification)))%if different numbers of channels differed across recordings in session, error out
                warning('Rercordings contain different Amplifications - per .meta files.  Will use number from first recording.')
            end
            DatsMetadata.Parameters.Amplification = Amplification(1);

            DatsMetadata.Parameters.VoltsPerUnit = VoltsPerUnit_Amplirec(basename,basepath);%calculate from .metas/.inis
            DatsMetadata.Parameters.BitsPerSample = 16;%amplirec default
            DatsMetadata.Parameters.VoltageRange = 10;%not used except to make xml
        end
    case 'Intan' %if not amplipex, look for intan
        d = drhd;
        NAmpChannels = [];
        NAuxChannels = [];
        NDigitalChannels = [];
        SamplingRate = [];
        thisidx = 0; 
        
        for didx = 1:length(d)
            if d(didx).isdir
                t = dir(fullfile(basepath,d(didx).name,'info.rhd'));
                if ~isempty(t)%only use folders with info.rhd inside

                   DatsMetadata.RecordingSystem = 'Intan';
                   rhd = fullfile(basepath,d(didx).name,'info.rhd');
    %                ampdats{end+1} = fullfile(basepath,d(didx).name,'amplifier.dat');
    %                recordingdiridxs = didx;
    
                   thisidx = thisidx + 1;
                  
                   DatsMetadata.Recordings.Names{thisidx} = d(didx).name;
                   
                   % read the info.rhd file for each individual recorded
                   % file/folder
                   try
                      [amplifier_channels, notes, aux_input_channels, spike_triggers,...         
                          board_dig_in_channels, supply_voltage_channels, frequency_parameters ] =...
                          read_Intan_RHD2000_file_In(fullfile(basepath,d(didx).name),'info.rhd');%function from Intan
                      DatsMetadata.Recordings.IntanRHDInfoRaw{thisidx} = v2struct(amplifier_channels, notes, aux_input_channels, spike_triggers,...         
                          board_dig_in_channels, supply_voltage_channels, frequency_parameters);
                      NAmpChannels(thisidx) = length(amplifier_channels);
                      NAuxChannels(thisidx) = length(aux_input_channels);
                      NDigitalChannels(thisidx) = length(board_dig_in_channels);
                      SamplingRate(thisidx) = frequency_parameters.amplifier_sample_rate;

                   catch%if read error, set values based on prior
                       if ~isempty(NAmpChannels)
%                            if ~isnan(NAmpChannels(1));
                          NAmpChannels(thisidx) = NAmpChannels(1);
                       else
                          sess = bz_getSessionInfo(basepath,'noPrompts',true);
                          NAmpChannels(thisidx) = sess.nChannels;
                       end
                       if ~isempty(NAuxChannels)
                           NAuxChannels(thisidx) = NAuxChannels(1);
                       else
                           NAuxChannels(thisidx) = nan;
                       end
                       if ~isempty(NDigitalChannels)
                           NDigitalChannels(thisidx) = NDigitalChannels(1);
                       else
                           NDigitalChannels(thisidx) = nan;
                       end
                       if ~isempty(SamplingRate)
                           SamplingRate(thisidx) = SamplingRate(1);
                       else
                           SamplingRate(thisidx) = nan;
                       end
                   end

                   t = dir(fullfile(basepath,d(didx).name,'amplifier.dat'));
                   t2 = dir(fullfile(basepath,d(didx).name,'amplifier_analogin_auxiliary_int16.dat'));
                   if ~isempty(t2) %if original .dat still around.  Handling Luke's RHD software to put all 16bit into same file
                       DatsMetadata.Recordings.Bytes(thisidx) = t2.bytes;                       
                   elseif ~isempty(t) %if original .dat still around
                       DatsMetadata.Recordings.Bytes(thisidx) = t.bytes;
                   else%if no .dat around
                       disp([DatsMetadata.Recordings.Names{end} ' not found']);
                       timedat = dir(fullfile(basepath,d(didx).name,'time.dat'));
                       if ~isnan(NAmpChannels(end));
                          DatsMetadata.Recordings.Bytes(thisidx) = timedat(1).bytes/2*NAmpChannels(end);
                       else
                           sess = bz_getSessionInfo(basepath);
                           DatsMetadata.Recordings.Bytes(thisidx) = timedat(1).bytes/2*sess.nChannels;
                       end
                   end

               end
            end
        end
        %handle errors created by improper read of rhd.info files
        if isnan(sum(abs(diff(NAmpChannels))))
            warning('At least one read error on rhd.info files.  Using values from other files to populate.')
            ok = find(~isnan(NAmpChannels),1,'first');
            if isempty(ok)
                warning('No good reads on NumberofAmpChannels. Assuming based on probe maps')
                sess = bz_getSessionInfo(basepath);
                NAmpChannels(isnan(NAmpChannels)) = sess.nChannels;
            else
                NAmpChannels(isnan(NAmpChannels)) = NAmpChannels(ok);
            end
        end
        if isnan(sum(abs(diff(NAuxChannels))))
            warning('At least one read error on rhd.info files.  Using values from other files to populate.')
            ok = find(~isnan(NAuxChannels),1,'first');
            if isempty(ok)
                error('No good reads on NAuxChannels, assuming 0')
                NAuxChannels(isnan(NAuxChannels)) = 0;
            else
                NAuxChannels(isnan(NAuxChannels)) = NAuxChannels(ok);
            end
        end
        if isnan(sum(abs(diff(NDigitalChannels))))
            warning('At least one read error on rhd.info files.  Using values from other files to populate.')
            ok = find(~isnan(NDigitalChannels),1,'first');
            if isempty(ok)
                error('No good reads on NDigitalChannels, assuming 0')
                NDigitalChannels(isnan(NDigitalChannels)) = 0;
            else
                NDigitalChannels(isnan(NDigitalChannels)) = NDigitalChannels(ok);
            end
        end
        if isnan(sum(abs(diff(SamplingRate))))
            warning('At least one read error on rhd.info files.  Using values from other files to populate.')
            ok = find(~isnan(SamplingRate),1,'first');
            if isempty(ok)
                error('No good reads on SamplingRate, taking from AnimalMetadata')
                sess = bz_getSessionInfo(basepath);
                SamplingRate(isnan(SamplingRate)) = sess.rates.wideband;
            else
                SamplingRate(isnan(SamplingRate)) = SamplingRate(ok);
            end
        end

        %handle errors where different parameters read for different files (rather than uniform)            
        if sum(abs(diff(NAmpChannels)))%if different numbers of channels differed across recordings in session, error out
            warning('Rercordings contain different numbers of channels - per info.rhd files.  Will use count from first recording.')
        end
        DatsMetadata.Parameters.NumberOfChannels = NAmpChannels(1);

        if sum(abs(diff(NAuxChannels)))%if different numbers of channels differed across recordings in session, error out
            warning('Rercordings contain different numbers auxiliary of channels - per info.rhd files.  Will use count from first recording.')
        end
        DatsMetadata.NumAuxiliaryChannels = NAuxChannels(1);

        if sum(abs(diff(NDigitalChannels)))%if different numbers of channels differed across recordings in session, error out
            warning('Rercordings contain different numbers digital of channels - per info.rhd files.  Will use count from first recording.')
        end
        DatsMetadata.NumDigitalChannels = NAuxChannels(1);

        %sampling rate
        if sum(abs(diff(SamplingRate)))%if different numbers of channels differed across recordings in session, error out
            warning('Rercordings contain different SamplingRates - per info.rhd files.  Will use number from first recording.')
        end
        DatsMetadata.Parameters.SampleRate = SamplingRate(1);

        % Not gather-able from intan info.rhd files
        DatsMetadata.Parameters.Amplification = 1;%digitized on chip, let's say 1
        DatsMetadata.Parameters.VoltsPerUnit = 0.0000002;%intan default                
        DatsMetadata.Parameters.BitsPerSample = 16;%intan default
        DatsMetadata.Parameters.VoltageRange = 10;%not used except to make xml

        %Now that everything is in order, calculate seconds from bytes
        DatsMetadata.Recordings.Seconds = DatsMetadata.Recordings.Bytes/DatsMetadata.Parameters.NumberOfChannels/2/DatsMetadata.Parameters.SampleRate;        
    case '' %if no known recording system used
        warning('No recording system found, using defaults.  May manually set values in guts of this .m file')
        DatsMetadata.RecordingSystem = 'Unknown';
        DatsMetadata.Recordings.Names = {};
        DatsMetadata.Recordings.Bytes = [];
        DatsMetadata.Recordings.Seconds = [];%can manually write in here if desired

        try 
            sess = bz_getSessionInfo(basename);
        catch
            disp('No metadata or session info detected, quitting bz_DatFileMetadata')
            return
        end
        
        DatsMetadata.Parameters.NumberOfChannels = sess.nChannels;
        DatsMetadata.Parameters.SampleRate = sess.rates.wideband;%Lab default
        DatsMetadata.Parameters.Amplification = 1;%Intan digitized on chip, let's say 1
        DatsMetadata.Parameters.VoltsPerUnit = 0.0000002;%Intan default                
        DatsMetadata.Parameters.BitsPerSample = sess.nBits;%Intan default
        DatsMetadata.Parameters.VoltageRange = 10;%not used except to make xml
end
    
save(fullfile(basepath,[basename '_DatsMetadata.mat']),'DatsMetadata')


%%
function [amplifier_channels, notes, aux_input_channels, spike_triggers,...         
board_dig_in_channels, supply_voltage_channels, frequency_parameters ] =...
read_Intan_RHD2000_file_In(path,file)

% read_Intan_RHD2000_file
%
% Version 1.3, 10 December 2013
%
% Reads Intan Technologies RHD2000 data file generated by evaluation board
% GUI.  Data are parsed and placed into variables that appear in the base
% MATLAB workspace.  Therefore, it is recommended to execute a 'clear'
% command before running this program to clear all other variables from the
% base workspace.
%
% INPUTS 
%   path and file are the path and filename for the rhd to be read.
%
% Example:
% >> clear
% >> read_Intan_RHD200_file
% >> whos
% >> amplifier_channels(1)
% >> plot(t_amplifier, amplifier_data(1,:))

if ~exist('path','var')
    [file, path, filterindex] = ...
        uigetfile('*.rhd', 'Select an RHD2000 Data File', 'MultiSelect', 'off');
else
    if ~strcmp('.rhd',file(end-3:end))
        file = strcat(file,'.rhd');
    end
end

% Read most recent file automatically.
%path = 'C:\Users\Reid\Documents\RHD2132\testing\';
%d = dir([path '*.rhd']);
%file = d(end).name;

tic;
filename = fullfile(path,file);
fid = fopen(filename, 'r');

s = dir(filename);
filesize = s.bytes;

% Check 'magic number' at beginning of file to make sure this is an Intan
% Technologies RHD2000 data file.
magic_number = fread(fid, 1, 'uint32');
if magic_number ~= hex2dec('c6912702')
    error('Unrecognized file type.');
end

% Read version number.
data_file_main_version_number = fread(fid, 1, 'int16');
data_file_secondary_version_number = fread(fid, 1, 'int16');

% fprintf(1, '\n');
% fprintf(1, 'Reading Intan Technologies RHD2000 Data File, Version %d.%d\n', ...
%     data_file_main_version_number, data_file_secondary_version_number);
% fprintf(1, '\n');

% Read information of sampling rate and amplifier frequency settings.
sample_rate = fread(fid, 1, 'single');
dsp_enabled = fread(fid, 1, 'int16');
actual_dsp_cutoff_frequency = fread(fid, 1, 'single');
actual_lower_bandwidth = fread(fid, 1, 'single');
actual_upper_bandwidth = fread(fid, 1, 'single');

desired_dsp_cutoff_frequency = fread(fid, 1, 'single');
desired_lower_bandwidth = fread(fid, 1, 'single');
desired_upper_bandwidth = fread(fid, 1, 'single');

% This tells us if a software 50/60 Hz notch filter was enabled during
% the data acquisition.
notch_filter_mode = fread(fid, 1, 'int16');
notch_filter_frequency = 0;
if (notch_filter_mode == 1)
    notch_filter_frequency = 50;
elseif (notch_filter_mode == 2)
    notch_filter_frequency = 60;
end

desired_impedance_test_frequency = fread(fid, 1, 'single');
actual_impedance_test_frequency = fread(fid, 1, 'single');

% Place notes in data strucure
notes = struct( ...
    'note1', fread_QString_In(fid), ...
    'note2', fread_QString_In(fid), ...
    'note3', fread_QString_In(fid) );
    
% If data file is from GUI v1.1 or later, see if temperature sensor data
% was saved.
num_temp_sensor_channels = 0;
if ((data_file_main_version_number == 1 & data_file_secondary_version_number >= 1) ...
    | (data_file_main_version_number > 1))
    num_temp_sensor_channels = fread(fid, 1, 'int16');
end

% If data file is from GUI v1.3 or later, load eval board mode.
eval_board_mode = 0;
if ((data_file_main_version_number == 1 & data_file_secondary_version_number >= 3) ...
    | (data_file_main_version_number > 1))
    eval_board_mode = fread(fid, 1, 'int16');
end

% Place frequency-related information in data structure.
frequency_parameters = struct( ...
    'amplifier_sample_rate', sample_rate, ...
    'aux_input_sample_rate', sample_rate / 4, ...
    'supply_voltage_sample_rate', sample_rate / 60, ...
    'board_adc_sample_rate', sample_rate, ...
    'board_dig_in_sample_rate', sample_rate, ...
    'desired_dsp_cutoff_frequency', desired_dsp_cutoff_frequency, ...
    'actual_dsp_cutoff_frequency', actual_dsp_cutoff_frequency, ...
    'dsp_enabled', dsp_enabled, ...
    'desired_lower_bandwidth', desired_lower_bandwidth, ...
    'actual_lower_bandwidth', actual_lower_bandwidth, ...
    'desired_upper_bandwidth', desired_upper_bandwidth, ...
    'actual_upper_bandwidth', actual_upper_bandwidth, ...
    'notch_filter_frequency', notch_filter_frequency, ...
    'desired_impedance_test_frequency', desired_impedance_test_frequency, ...
    'actual_impedance_test_frequency', actual_impedance_test_frequency );

% Define data structure for spike trigger settings.
spike_trigger_struct = struct( ...
    'voltage_trigger_mode', {}, ...
    'voltage_threshold', {}, ...
    'digital_trigger_channel', {}, ...
    'digital_edge_polarity', {} );

new_trigger_channel = struct(spike_trigger_struct);
spike_triggers = struct(spike_trigger_struct);

% Define data structure for data channels.
channel_struct = struct( ...
    'native_channel_name', {}, ...
    'custom_channel_name', {}, ...
    'native_order', {}, ...
    'custom_order', {}, ...
    'board_stream', {}, ...
    'chip_channel', {}, ...
    'port_name', {}, ...
    'port_prefix', {}, ...
    'port_number', {}, ...
    'electrode_impedance_magnitude', {}, ...
    'electrode_impedance_phase', {} );

new_channel = struct(channel_struct);

% Create structure arrays for each type of data channel.
amplifier_channels = struct(channel_struct);
aux_input_channels = struct(channel_struct);
supply_voltage_channels = struct(channel_struct);
board_adc_channels = struct(channel_struct);
board_dig_in_channels = struct(channel_struct);
board_dig_out_channels = struct(channel_struct);

amplifier_index = 1;
aux_input_index = 1;
supply_voltage_index = 1;
board_adc_index = 1;
board_dig_in_index = 1;
board_dig_out_index = 1;

% Read signal summary from data file header.

number_of_signal_groups = fread(fid, 1, 'int16');

for signal_group = 1:number_of_signal_groups
    signal_group_name = fread_QString_In(fid);
    signal_group_prefix = fread_QString_In(fid);
    signal_group_enabled = fread(fid, 1, 'int16');
    signal_group_num_channels = fread(fid, 1, 'int16');
    signal_group_num_amp_channels = fread(fid, 1, 'int16');

    if (signal_group_num_channels > 0 && signal_group_enabled > 0)
        new_channel(1).port_name = signal_group_name;
        new_channel(1).port_prefix = signal_group_prefix;
        new_channel(1).port_number = signal_group;
        for signal_channel = 1:signal_group_num_channels
            new_channel(1).native_channel_name = fread_QString_In(fid);
            new_channel(1).custom_channel_name = fread_QString_In(fid);
            new_channel(1).native_order = fread(fid, 1, 'int16');
            new_channel(1).custom_order = fread(fid, 1, 'int16');
            signal_type = fread(fid, 1, 'int16');
            channel_enabled = fread(fid, 1, 'int16');
            new_channel(1).chip_channel = fread(fid, 1, 'int16');
            new_channel(1).board_stream = fread(fid, 1, 'int16');
            new_trigger_channel(1).voltage_trigger_mode = fread(fid, 1, 'int16');
            new_trigger_channel(1).voltage_threshold = fread(fid, 1, 'int16');
            new_trigger_channel(1).digital_trigger_channel = fread(fid, 1, 'int16');
            new_trigger_channel(1).digital_edge_polarity = fread(fid, 1, 'int16');
            new_channel(1).electrode_impedance_magnitude = fread(fid, 1, 'single');
            new_channel(1).electrode_impedance_phase = fread(fid, 1, 'single');
            
            if (channel_enabled)
                switch (signal_type)
                    case 0
                        amplifier_channels(amplifier_index) = new_channel;
                        spike_triggers(amplifier_index) = new_trigger_channel;
                        amplifier_index = amplifier_index + 1;
                    case 1
                        aux_input_channels(aux_input_index) = new_channel;
                        aux_input_index = aux_input_index + 1;
                    case 2
                        supply_voltage_channels(supply_voltage_index) = new_channel;
                        supply_voltage_index = supply_voltage_index + 1;
                    case 3
                        board_adc_channels(board_adc_index) = new_channel;
                        board_adc_index = board_adc_index + 1;
                    case 4
                        board_dig_in_channels(board_dig_in_index) = new_channel;
                        board_dig_in_index = board_dig_in_index + 1;
                    case 5
                        board_dig_out_channels(board_dig_out_index) = new_channel;
                        board_dig_out_index = board_dig_out_index + 1;
                    otherwise
                        error('Unknown channel type');
                end
            end
            
        end
    end
end

% Summarize contents of data file.
num_amplifier_channels = amplifier_index - 1;
num_aux_input_channels = aux_input_index - 1;
num_supply_voltage_channels = supply_voltage_index - 1;
num_board_adc_channels = board_adc_index - 1;
num_board_dig_in_channels = board_dig_in_index - 1;
num_board_dig_out_channels = board_dig_out_index - 1;

% fprintf(1, 'Found %d amplifier channel%s.\n', ...
%     num_amplifier_channels, plural(num_amplifier_channels));
% fprintf(1, 'Found %d auxiliary input channel%s.\n', ...
%     num_aux_input_channels, plural(num_aux_input_channels));
% fprintf(1, 'Found %d supply voltage channel%s.\n', ...
%     num_supply_voltage_channels, plural(num_supply_voltage_channels));
% fprintf(1, 'Found %d board ADC channel%s.\n', ...
%     num_board_adc_channels, plural(num_board_adc_channels));
% fprintf(1, 'Found %d board digital input channel%s.\n', ...
%     num_board_dig_in_channels, plural(num_board_dig_in_channels));
% fprintf(1, 'Found %d board digital output channel%s.\n', ...
%     num_board_dig_out_channels, plural(num_board_dig_out_channels));
% fprintf(1, 'Found %d temperature sensors channel%s.\n', ...
%     num_temp_sensor_channels, plural(num_temp_sensor_channels));
% fprintf(1, '\n');

% Determine how many samples the data file contains.

% Each data block contains 60 amplifier samples.
bytes_per_block = 60 * 4;  % timestamp data
bytes_per_block = bytes_per_block + 60 * 2 * num_amplifier_channels;
% Auxiliary inputs are sampled 4x slower than amplifiers
bytes_per_block = bytes_per_block + 15 * 2 * num_aux_input_channels;
% Supply voltage is sampled 60x slower than amplifiers
bytes_per_block = bytes_per_block + 1 * 2 * num_supply_voltage_channels;
% Board analog inputs are sampled at same rate as amplifiers
bytes_per_block = bytes_per_block + 60 * 2 * num_board_adc_channels;
% Board digital inputs are sampled at same rate as amplifiers
if (num_board_dig_in_channels > 0)
    bytes_per_block = bytes_per_block + 60 * 2;
end
% Board digital outputs are sampled at same rate as amplifiers
if (num_board_dig_out_channels > 0)
    bytes_per_block = bytes_per_block + 60 * 2;
end
% Temp sensor is sampled 60x slower than amplifiers
if (num_temp_sensor_channels > 0)
   bytes_per_block = bytes_per_block + 1 * 2 * num_temp_sensor_channels; 
end

% How many data blocks remain in this file?
data_present = 0;
bytes_remaining = filesize - ftell(fid);
if (bytes_remaining > 0)
    data_present = 1;
end

num_data_blocks = bytes_remaining / bytes_per_block;

num_amplifier_samples = 60 * num_data_blocks;
num_aux_input_samples = 15 * num_data_blocks;
num_supply_voltage_samples = 1 * num_data_blocks;
num_board_adc_samples = 60 * num_data_blocks;
num_board_dig_in_samples = 60 * num_data_blocks;
num_board_dig_out_samples = 60 * num_data_blocks;

record_time = num_amplifier_samples / sample_rate;

% if (data_present)
%     fprintf(1, 'File contains %0.3f seconds of data.  Amplifiers were sampled at %0.2f kS/s.\n', ...
%         record_time, sample_rate / 1000);
%     fprintf(1, '\n');
% else
%     fprintf(1, 'Header file contains no data.  Amplifiers were sampled at %0.2f kS/s.\n', ...
%         sample_rate / 1000);
%     fprintf(1, '\n');
% end

if (data_present)
    
    % Pre-allocate memory for data.
    fprintf(1, 'Allocating memory for data...\n');

    t_amplifier = zeros(1, num_amplifier_samples);

    amplifier_data = zeros(num_amplifier_channels, num_amplifier_samples);
    aux_input_data = zeros(num_aux_input_channels, num_aux_input_samples);
    supply_voltage_data = zeros(num_supply_voltage_channels, num_supply_voltage_samples);
    temp_sensor_data = zeros(num_temp_sensor_channels, num_supply_voltage_samples);
    board_adc_data = zeros(num_board_adc_channels, num_board_adc_samples);
    board_dig_in_data = zeros(num_board_dig_in_channels, num_board_dig_in_samples);
    board_dig_in_raw = zeros(1, num_board_dig_in_samples);
    board_dig_out_data = zeros(num_board_dig_out_channels, num_board_dig_out_samples);
    board_dig_out_raw = zeros(1, num_board_dig_out_samples);

    % Read sampled data from file.
    fprintf(1, 'Reading data from file...\n');

    amplifier_index = 1;
    aux_input_index = 1;
    supply_voltage_index = 1;
    board_adc_index = 1;
    board_dig_in_index = 1;
    board_dig_out_index = 1;

    print_increment = 10;
    percent_done = print_increment;
    for i=1:num_data_blocks
        % In version 1.2, we moved from saving timestamps as unsigned
        % integeters to signed integers to accomidate negative (adjusted)
        % timestamps for pretrigger data.
        if ((data_file_main_version_number == 1 && data_file_secondary_version_number >= 2) ...
        || (data_file_main_version_number > 1))
            t_amplifier(amplifier_index:(amplifier_index+59)) = fread(fid, 60, 'int32');
        else
            t_amplifier(amplifier_index:(amplifier_index+59)) = fread(fid, 60, 'uint32');
        end
        if (num_amplifier_channels > 0)
            amplifier_data(:, amplifier_index:(amplifier_index+59)) = fread(fid, [60, num_amplifier_channels], 'uint16')';
        end
        if (num_aux_input_channels > 0)
            aux_input_data(:, aux_input_index:(aux_input_index+14)) = fread(fid, [15, num_aux_input_channels], 'uint16')';
        end
        if (num_supply_voltage_channels > 0)
            supply_voltage_data(:, supply_voltage_index) = fread(fid, [1, num_supply_voltage_channels], 'uint16')';
        end
        if (num_temp_sensor_channels > 0)
            temp_sensor_data(:, supply_voltage_index) = fread(fid, [1, num_temp_sensor_channels], 'int16')';
        end
        if (num_board_adc_channels > 0)
            board_adc_data(:, board_adc_index:(board_adc_index+59)) = fread(fid, [60, num_board_adc_channels], 'uint16')';
        end
        if (num_board_dig_in_channels > 0)
            board_dig_in_raw(board_dig_in_index:(board_dig_in_index+59)) = fread(fid, 60, 'uint16');
        end
        if (num_board_dig_out_channels > 0)
            board_dig_out_raw(board_dig_out_index:(board_dig_out_index+59)) = fread(fid, 60, 'uint16');
        end

        amplifier_index = amplifier_index + 60;
        aux_input_index = aux_input_index + 15;
        supply_voltage_index = supply_voltage_index + 1;
        board_adc_index = board_adc_index + 60;
        board_dig_in_index = board_dig_in_index + 60;
        board_dig_out_index = board_dig_out_index + 60;

        fraction_done = 100 * (i / num_data_blocks);
        if (fraction_done >= percent_done)
            fprintf(1, '%d%% done...\n', percent_done);
            percent_done = percent_done + print_increment;
        end
    end

    % Make sure we have read exactly the right amount of data.
    bytes_remaining = filesize - ftell(fid);
    if (bytes_remaining ~= 0)
        %error('Error: End of file not reached.');
    end

end

% Close data file.
fclose(fid);

if (data_present)
    
    fprintf(1, 'Parsing data...\n');

    % Extract digital input channels to separate variables.
    for i=1:num_board_dig_in_channels
       mask = 2^(board_dig_in_channels(i).native_order) * ones(size(board_dig_in_raw));
       board_dig_in_data(i, :) = (bitand(board_dig_in_raw, mask) > 0);
    end
    for i=1:num_board_dig_out_channels
       mask = 2^(board_dig_out_channels(i).native_order) * ones(size(board_dig_out_raw));
       board_dig_out_data(i, :) = (bitand(board_dig_out_raw, mask) > 0);
    end

    % Scale voltage levels appropriately.
    amplifier_data = 0.195 * (amplifier_data - 32768); % units = microvolts
    aux_input_data = 37.4e-6 * aux_input_data; % units = volts
    supply_voltage_data = 74.8e-6 * supply_voltage_data; % units = volts
    if (eval_board_mode == 1)
        board_adc_data = 152.59e-6 * (board_adc_data - 32768); % units = volts    
    else
        board_adc_data = 50.354e-6 * board_adc_data; % units = volts
    end
    temp_sensor_data = temp_sensor_data / 100; % units = deg C

    % Check for gaps in timestamps.
    num_gaps = sum(diff(t_amplifier) ~= 1);
    if (num_gaps == 0)
        fprintf(1, 'No missing timestamps in data.\n');
    else
        fprintf(1, 'Warning: %d gaps in timestamp data found.  Time scale will not be uniform!\n', ...
            num_gaps);
    end

    % Scale time steps (units = seconds).
    t_amplifier = t_amplifier / sample_rate;
    t_aux_input = t_amplifier(1:4:end);
    t_supply_voltage = t_amplifier(1:60:end);
    t_board_adc = t_amplifier;
    t_dig = t_amplifier;
    t_temp_sensor = t_supply_voltage;

    % If the software notch filter was selected during the recording, apply the
    % same notch filter to amplifier data here.
    if (notch_filter_frequency > 0)
        fprintf(1, 'Applying notch filter...\n');

        print_increment = 10;
        percent_done = print_increment;
        for i=1:num_amplifier_channels
            amplifier_data(i,:) = ...
                notch_filter_In(amplifier_data(i,:), sample_rate, notch_filter_frequency, 10);

            fraction_done = 100 * (i / num_amplifier_channels);
            if (fraction_done >= percent_done)
                fprintf(1, '%d%% done...\n', percent_done);
                percent_done = percent_done + print_increment;
            end

        end
    end

end

% Move variables to base workspace.

move_to_base_workspace_In(notes);
move_to_base_workspace_In(frequency_parameters);

if (num_amplifier_channels > 0)
    move_to_base_workspace_In(amplifier_channels);
    if (data_present)
        move_to_base_workspace_In(amplifier_data);
        move_to_base_workspace_In(t_amplifier);
    end
    move_to_base_workspace_In(spike_triggers);
end
if (num_aux_input_channels > 0)
    move_to_base_workspace_In(aux_input_channels);
    if (data_present)
        move_to_base_workspace_In(aux_input_data);
        move_to_base_workspace_In(t_aux_input);
    end
end
if (num_supply_voltage_channels > 0)
    move_to_base_workspace_In(supply_voltage_channels);
    if (data_present)
        move_to_base_workspace_In(supply_voltage_data);
        move_to_base_workspace_In(t_supply_voltage);
    end
end
if (num_board_adc_channels > 0)
    move_to_base_workspace_In(board_adc_channels);
    if (data_present)
        move_to_base_workspace_In(board_adc_data);
        move_to_base_workspace_In(t_board_adc);
    end
end
if (num_board_dig_in_channels > 0)
    move_to_base_workspace_In(board_dig_in_channels);
    if (data_present)
        move_to_base_workspace_In(board_dig_in_data);
        move_to_base_workspace_In(t_dig);
    end
end
if (num_board_dig_out_channels > 0)
    move_to_base_workspace_In(board_dig_out_channels);
    if (data_present)
        move_to_base_workspace_In(board_dig_out_data);
        move_to_base_workspace_In(t_dig);
    end
end
if (num_temp_sensor_channels > 0)
    if (data_present)
        move_to_base_workspace_In(temp_sensor_data);
        move_to_base_workspace_In(t_temp_sensor);
    end
end

% fprintf(1, 'Done!  Elapsed time: %0.1f seconds\n', toc);
% if (data_present)
%     fprintf(1, 'Extracted data are now available in the MATLAB workspace.\n');
% else
%     fprintf(1, 'Extracted waveform information is now available in the MATLAB workspace.\n');
% end
% fprintf(1, 'Type ''whos'' to see variables.\n');
% fprintf(1, '\n');

return


function a = fread_QString_In(fid)

% a = read_QString(fid)
%
% Read Qt style QString.  The first 32-bit unsigned number indicates
% the length of the string (in bytes).  If this number equals 0xFFFFFFFF,
% the string is null.

a = '';
length = fread(fid, 1, 'uint32');
if length == hex2num('ffffffff')
    return;
end
% convert length from bytes to 16-bit Unicode words
length = length / 2;

for i=1:length
    a(i) = fread(fid, 1, 'uint16');
end

return


function s = plural_In(n)

% s = plural(n)
% 
% Utility function to optionally plurailze words based on the value
% of n.

if (n == 1)
    s = '';
else
    s = 's';
end

return


function out = notch_filter_In(in, fSample, fNotch, Bandwidth)

% out = notch_filter(in, fSample, fNotch, Bandwidth)
%
% Implements a notch filter (e.g., for 50 or 60 Hz) on vector 'in'.
% fSample = sample rate of data (in Hz or Samples/sec)
% fNotch = filter notch frequency (in Hz)
% Bandwidth = notch 3-dB bandwidth (in Hz).  A bandwidth of 10 Hz is
%   recommended for 50 or 60 Hz notch filters; narrower bandwidths lead to
%   poor time-domain properties with an extended ringing response to
%   transient disturbances.
%
% Example:  If neural data was sampled at 30 kSamples/sec
% and you wish to implement a 60 Hz notch filter:
%
% out = notch_filter(in, 30000, 60, 10);

tstep = 1/fSample;
Fc = fNotch*tstep;

L = length(in);

% Calculate IIR filter parameters
d = exp(-2*pi*(Bandwidth/2)*tstep);
b = (1 + d*d)*cos(2*pi*Fc);
a0 = 1;
a1 = -b;
a2 = d*d;
a = (1 + d*d)/2;
b0 = 1;
b1 = -2*cos(2*pi*Fc);
b2 = 1;

out = zeros(size(in));
out(1) = in(1);  
out(2) = in(2);
% (If filtering a continuous data stream, change out(1) and out(2) to the
%  previous final two values of out.)

% Run filter
for i=3:L
    out(i) = (a*b2*in(i-2) + a*b1*in(i-1) + a*b0*in(i) - a2*out(i-2) - a1*out(i-1))/a0;
end

return


function move_to_base_workspace_In(variable)

% move_to_base_workspace(variable)
%
% Move variable from function workspace to base MATLAB workspace so
% user will have access to it after the program ends.

variable_name = inputname(1);
assignin('base', variable_name, variable);

return;