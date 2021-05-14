function accel = bz_getIntanAccel(varargin)
% bz_getIntanAccel - Get accelerometer data from Intan-generated auxiliary.dat file
%
% USAGE
%
%    acc = bz_getIntanAccel(varargin)
%
% INPUTS
%
%    basepath        -path to recording (where .dat/.clu/etc files are)
%    lowpass         -numeric value to low-pass filter acceleration data           
%    samplingRate    -numeric value to control final sampling rate of acceleration data            
%    forceReload     -logical (default=false) to force loading from auxiliary.dat file            
%    saveMat         -logical (default=false) to save in buzcode format
%    noPrompts       -logical (default=false) to supress any user prompts
%
% OUTPUTS
%
%    acc - behavior struct with the following fields
%          .timestamps     -array of timestamps that match the data subfields (in seconds)
%          .acceleration   -data substruct with x, y, z acceleration and overall magnitude of acceleration
%          .samplingRate   -sampling rate of data
%          .lowpass        -value used to low-pass filter raw data   
%          .units          -unit of measurement of acceleration (volts, see the note below)
%          .behaviorinfo   -information substruct
% 
% NOTES:
% 
% Extract accelerometer data from auxiliary.dat file - when doing this, downsample 
% to the rate of LFP acquisition, and low-pass filter the resulting signals 
% (parameter 'lowpass'). Use acceleration along x, y, and z axes to compute 
% magnitude of acceleration (i.e., the norm of a 3D vector). This can be used 
% as a proxy for whether the animal  is moving or not. In a final step, average 
% subintervals of acceleration values to get a desired number of values per 
% second (parameter 'samplingRate'). Use info.rhd Intan data file to verify
% number of active auxiliary.dat files - this step depends on a modified
% version of 'read_Intan_RHD2000_file.m' (Intan MATLAB script)
%
% The Intan accelerometer stores voltage values that reflect instantaneous
% acceleration along the x, y, and z axes. On average, a value of 0.34 V maps
% to 1g acceleration. For details, see the following link:
% http://intantech.com/files/Intan_RHD2000_accelerometer_calibration.pdf
% 
% Written by Roman Huszar, 2019

%% Process and check user input

% Parse user input
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'lowpass', 1, @isnumeric);
addParameter(p,'samplingRate', 1, @isnumeric);
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'noPrompts',false,@islogical);
parse(p,varargin{:})

% Store user input
basepath = p.Results.basepath;
lowpass = p.Results.lowpass;
samplingRate = p.Results.samplingRate;
forceReload = p.Results.forceReload;
saveMat = p.Results.saveMat;
noPrompts = p.Results.noPrompts;

% Session information
basename = bz_BasenameFromBasepath(basepath);
aux_input_path = fullfile(basepath, [basename '_auxiliary.dat']);
accel_output_path = fullfile(basepath, [basename '.accel.behavior.mat']);
sessionInfo = bz_getSessionInfo(basepath, 'noPrompts', noPrompts);
fs_wide = sessionInfo.rates.wideband;
fs_lfp  = sessionInfo.lfpSampleRate;

% If auxiliary file not find, crash with grace
if 2 ~= exist(aux_input_path, 'file')
    error('auxiliary.dat file not found in basepath - cannot extract accelerometer')
end

% If file exists, return its contents
if 2 == exist(accel_output_path, 'file') && ~forceReload
    disp('Loading acceleration data from acceler.behavior  file..')
    load(accel_output_path)   %#ok<LOAD>
    return
end

% Ask user about saving data
if ~noPrompts && ~saveMat 
    savebutton = questdlg(['Would you like to store acceleration ',...
                            'data in .accel.behavior.mat? ']);
    if strcmp(savebutton,'Yes'); saveMat = true; end
end

% Check validity of target sampling rate
if samplingRate > fs_lfp
    error('samplingRate must be smaller than lfp sampling rate, which is %d.\nType ''help bz_getIntanAccel'' for more.', fs_lfp)
end

% Get the auxiliary channels that are active
basepath_contents = dir(basepath);
dirs_in_basepath = {basepath_contents([basepath_contents.isdir]).name};
% Find the first 'info.rhd' file. This assumes that sessions with more dat files 
% (i.e., due to Intan crash) have identical info.rhd parameter settings...
% IMPORTANT NOTE - we assume info.rhd lives either in basepath directory,
% or in one of its subdirectories (we don't check sub-subdirectories!!)
found_flag = false;
for ii = 1:length(dirs_in_basepath)     % First two directories are '.' and '..'
    inforhd_path = fullfile(basepath, dirs_in_basepath{ii});
    if 2 == exist(fullfile(inforhd_path,'info.rhd'), 'file'); found_flag = true; break; end   % Exit when first info.rhd file found
end

% Read the Intan info file
if found_flag
    try
        intaninfo = read_Intan_RHD2000_file('returnStruct', true, 'filepath', inforhd_path);
        n_active_channels = length(intaninfo.aux_input_channels);
        active_channel_ids = {intaninfo.aux_input_channels.native_channel_name};
        is_x_on = any( cellfun(@(x) ~isempty(x), regexp(active_channel_ids, 'AUX1')) );
        is_y_on = any( cellfun(@(x) ~isempty(x), regexp(active_channel_ids, 'AUX2')) );
        is_z_on = any( cellfun(@(x) ~isempty(x), regexp(active_channel_ids, 'AUX3')) );
    catch
        warning('found info.rhd, but it is unreadable. Cannot confirm number of active channels in auxiliary.dat. Default is set to 3 - can cause problems!')
        n_active_channels = 3;
        is_x_on = true; is_y_on = true; is_z_on = true;
    end
else
    warning('found no info.rhd file in basepath subdirectories. Cannot confirm number of active channels in auxiliary.dat. Default is set to 3 - can cause problems!')
    n_active_channels = 3;
    is_x_on = true; is_y_on = true; is_z_on = true;
end
    
% Get IDs of the specific active auxiliary channels
if n_active_channels < 3
    warning('auxiliary.dat file contains fewer than the typical 3 channels, motion detection may be less accurate');
end

%% Process auxiliary dat file

% Load auxiliary data file - downsample to rate of LFP while doing so
disp('Loading auxiliary.dat binary file...')
aux_data = bz_LoadBinary(aux_input_path, 'nchannels', n_active_channels, 'channels', 1:n_active_channels, 'frequency', fs_wide, 'precision', 'uint16', 'downsample', round(fs_wide / fs_lfp));
aux_timestamps = [0 : 1/fs_lfp : (length(aux_data)-1)/fs_lfp]';

% Compute acceleration vector 
% (Take into account older versions of MATLAB)
try
    meeg = vecnorm(double(aux_data),2,2);
catch
    meeg = vecnorm(double(aux_data),2);
end
    
% Low pass filter the signal (default = 1 Hz) 
% NOTE: this assumes that finer timescale data are not needed  
Wn_theta = [0.1/(fs_lfp/2) lowpass/(fs_lfp/2)];
[btheta, atheta] = butter(2, Wn_theta);
meeg_filt = filtfilt(btheta, atheta, meeg);
meeg_filt = abs(meeg_filt);
aux_data_filt = filtfilt(btheta, atheta, double(aux_data));

% Average values in prespecified interval to get desired number of data
% points per second (samplingRate)
dat_per_samp = round(fs_lfp / samplingRate);
% Motion proxy
motion = mean(reshape(meeg_filt(1:(length(meeg_filt) - mod(length(meeg_filt), dat_per_samp))), dat_per_samp, []), 1);
% Each axis separately
ii = 1;
if is_x_on
    x = mean(reshape(aux_data_filt(1:(length(aux_data_filt) - mod(length(aux_data_filt), dat_per_samp)), ii), dat_per_samp, []), 1); ii = ii + 1;
else
    x = [];
end
if is_y_on
    y = mean(reshape(aux_data_filt(1:(length(aux_data_filt) - mod(length(aux_data_filt), dat_per_samp)), ii), dat_per_samp, []), 1); ii = ii + 1;
else
    y = [];
end
if is_z_on
    z = mean(reshape(aux_data_filt(1:(length(aux_data_filt) - mod(length(aux_data_filt), dat_per_samp)), 3), dat_per_samp, []), 1);
else
    z = [];
end

% Generate a buzcode behavior struct to hold the accelerometer information
accel = struct();
accel.timestamps                      = aux_timestamps(1:dat_per_samp:end-dat_per_samp);
accel.acceleration.x                  = x;
accel.acceleration.y                  = y;
accel.acceleration.z                  = z;
accel.acceleration.motion             = motion;
accel.acceleration.activeChannels     = n_active_channels;
accel.samplingRate                    = samplingRate;
accel.lowpass                         = lowpass;
accel.units                           = 'V';
accel.behaviorinfo.description        = 'accelerometer from Intan auxiliary dat file';
accel.behaviorinfo.acquisitionsystem  = 'Intan';
accel.behaviorinfo.processingfunction = 'bz_getIntanAccel.m';

% Save struct if prompted by user
if saveMat
    save(accel_output_path, 'accel')
end


end




