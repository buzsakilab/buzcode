function saved_binned_data_file_name = create_binned_data_from_raster_data(raster_file_directory_name, save_prefix_name, bin_width, sampling_interval, start_time, end_time)

% This function takes the name of a directory that contains files in raster-format 
%  and creates data that is in binned-format data that contain the data averaged
%  over some time period specified in the variable bin_wdith.  The arguments to this function are:
%
%  1. raster_file_directory_name: the path to the directory that contains the files
%      in raster format. The end of the path name can also contain a character string and only 
%      raster files that start with this character string will be included in the binned_data 
%      (e.g., if you have my_directory/*PFC* only raster files that contain PFC will be included 
%      in the binned_data. 
%
%  2. save_prefix_name: the beginning of a file name (possibly including a directory name)
%      that specifies the name that the binned data should be saved as.  Appended on 
%      to the end of this name when the file is saved is the bin width, step size, 
%      and possibly start and end times used in the binning.my_raster_file_directory/ 
%      
%  3. bin_width:  the bin size that is averaged over when creating binned-format features
%      (e.g., if the raster file has spike times given with millisecond position, and
%      bin_width = 500, then the binned-format data will contain average firing rates
%      (i.e., the spike-count rate) in 500 ms bins).  If bin_width is a vector that is the same size 
%      as sampling_interval then different bin widths can be specified for different time periods 
%      (e.g., if bin_widths = [50 500, 100], then the first bin will average over 50 ms of data, 
%      the second over 500 ms of data, and the third over 100 ms of data).
%      
%  4. sampling_interval: specifies the sampling interval between successive binned-data points
%     (e.g., if the raster file has spike times given with millisecond position, 
%     and if sampling_interval = 50, then a binned data point will be computed at 50ms intervals.
%     If sampling_interval is a vector, then the values in this vector will specify the start times
%     for the bins (e.g., if sampling_interval = [1 200 500], then the first bin will start time 1,
%     the second bin will start at time 200, and the third bin will start at time 500).
%
%  Optional arguments:
%
%  5. start_time:  This specifies the time to start the binning process.  If 
%      this argument is not set, then the binning will start with the first 
%      data point in the raster-file.
%
%  6. end_time:  This specifies the time when to end the binning process. If 
%      this argument is not set, then the binning will end with the last  
%      data point in the raster-file.
%
%  Note: if sampling_interval is a vector, the arguments 5 and 6 cannot be set (this is
%      because the start bin times have already been exactly specified by the sampling_interval vector).
%
%  This function returns the name of the saved file in the varable 'saved_binned_data_file_name'. It also saves the 
%       binned-format data a file in the directory specified (the file is saved as a version 7.3 mat file). 
%       The saved file contains the three variables needed to conform to binned-format which are:  
%       'binned_data' which contains the data, 'binned_labels' which contains the labels, and 'binned_site_info'
%       which contains the any extra information about each site. The data in the three saved variables are extracted 
%       from the raster-format files that are contained in the input directory.  Additional information 
%       is saved in binned_site_info.binning_parameters that contains the bin_width, sampling interval, etc., that were used to create this binned data.
%
%
%  Example 1:  
%    
%   Suppose we had a directory called my_raster_file_directory/ that contained a number
%   of files in raster-data format from the spike times of neurons specified at 1 ms 
%   resolution.  Then running:
%
%   create_binned_data_from_raster_data('my_raster_file_directory/', 'my_save_dir/binned_data', 150, 50, 200, 1000);
%    
%   will create a file in the binned_data_150ms_bins_50ms_sampled_200start_time_1000end_time
%   that will be saved in the directory my_save_dir/.  This file will contain the average
%   firing rate calculated over 150 ms intervals for all neurons in the directory 
%   my_raster_file_directory/.  The firing rates will be sampled every 50 ms, starting
%   200 ms into the raster-file data and ending 1000 ms into the raster-file data.
%
%  Example 2:  
%    
%   Suppose we had a directory called my_raster_file_directory/ that contained a number
%   of files in raster-data format from the spike times of neurons specified at 1 ms 
%   resolution.  Then running:
%
%   create_binned_data_from_raster_data('my_raster_file_directory/', 'my_save_dir/binned_data_100ms_bins_plus_2_extra', [100 .* ones(10, 1) 250 250], [1:50:500 1 251]);
%    
%   will create a file in the binned_data_100ms_bins_plus_2_extra_custom_bins_custom_sampled
%   that will be saved in the directory my_save_dir/.  This file will contain the average
%   firing rate calculated over 100 ms intervals for all neurons in the directory 
%   my_raster_file_directory/ and will be sampled every 50 ms for the first 500 ms of data.  
%   Additionally, there will be two 250 ms bins at the end of the binned data that have the average 
%   firing rates over 1-250 ms into the trial and 251-500 ms into the trial. 


%==========================================================================

%     This code is part of the Neural Decoding Toolbox.
%     Copyright (C) 2011 by Ethan Meyers (emeyers@mit.edu)
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
%========================================================================== 




%  fix the directory name in case there it does not end with a slash
last_char = raster_file_directory_name(end);
if (strcmp(last_char, '/') + strcmp(last_char, '\') + strcmp(last_char, '*')) == 0
    raster_file_directory_name = [raster_file_directory_name '/']
end


% if the directory ends with a name a * then use only those file names that have the string up to the *
% for example my_raster_directory/*PFC* will use only the files that have PFC in the file name 
if strcmp(last_char, '*')
    raster_file_dir = dir(raster_file_directory_name);
    [t r] = strtok(raster_file_directory_name, '/');
    while ~isempty(r)
        [t r] = strtok(r, '/');
    end
    raster_file_directory_name = raster_file_directory_name(1:(end - length(t)));  % remove the end string from the directory name
else
    raster_file_dir = dir([raster_file_directory_name '*.mat']);
end



if isempty(raster_file_dir)
    if isempty(raster_file_dir)
        error('The directory name given does not exist');
    else
        error('The directory name given does not contain any .mat files');
    end
end


if nargin < 5
    start_time = 1;
end
if nargin < 6
    load([raster_file_directory_name raster_file_dir(1).name]);
    end_time = size(raster_data, 2);  
end


if (length(bin_width) == 1) && (length(sampling_interval) == 1)  % if a single bin width and step size have been specified, then create binned data that averaged data over bin_width sized bins, sampled at sampling_interval intervals
    the_bin_start_times = start_time:sampling_interval:(end_time - bin_width  + 1);
    the_bin_widths = bin_width .* ones(size(the_bin_start_times));  
elseif (length(bin_width) == 1) && (length(sampling_interval) > 1) % if multiple step size have been specified indicating bin start times, but only a single bin width is given, then create binned data that averaged data over bin_width sized bins, sampled using the start bin times given by sampling_interval
    the_bin_start_times = sampling_interval;
    the_bin_widths = bin_width .* ones(size(the_bin_start_times)); 
elseif (length(bin_width) > 1) && (length(sampling_interval) > 1)     
    if length(bin_width) ~= length(sampling_interval), error('If a series of bin widths are specified by having bin_width be a vector, and a series of start binning times are specified by having sampling_interval be a vector, then the length(bin_width) must equal length(sampling_interval)'); end
    the_bin_start_times = sampling_interval;
    the_bin_widths = bin_width;
end
    

if (length(sampling_interval) > 1) && (nargin > 4)
    error('If a series of start bin times are given by having sampling_interval be a vector, then one can not specify start_time and end_time arguments (i.e., the number of arguments to this function must be less than 5');
end


fprintf('\n');

% go through all the files and bin them
for i = 1:length(raster_file_dir)


    % print a message the the data is being binned (and add a dot for each file that has been binned
    curr_bin_string = [' Binning the data: ' num2str(i) ' of ' num2str(length(raster_file_dir))];
    if i == 1
        disp(curr_bin_string); 
    else
        fprintf([repmat(8,1,bin_str_len) curr_bin_string]);         
    end
    bin_str_len = length(curr_bin_string);

    
    
   load([raster_file_directory_name raster_file_dir(i).name]);  

   curr_binned_data = bin_one_site(raster_data, the_bin_start_times, the_bin_widths);  % use the below helper function to bin the data

   binned_data{i} = curr_binned_data;

   
   % save all the labels
   the_label_names = fieldnames(raster_labels);
   for iLabel = 1:length(the_label_names)
      eval(['binned_labels.' the_label_names{iLabel} '{i} = raster_labels.' the_label_names{iLabel} ';']); 
   end
       
   
   % save any extra neuron info 
   if ~isempty(raster_site_info)
       the_info_field_names = fieldnames(raster_site_info);  
       for iInfo = 1:length(the_info_field_names)
       
          if isstr( eval(['raster_site_info.' the_info_field_names{iInfo}])) 
              eval(['binned_site_info.' the_info_field_names{iInfo} '{i} = raster_site_info.' the_info_field_names{iInfo} ';']); 
          elseif isempty( eval(['raster_site_info.' the_info_field_names{iInfo}]))
             % just ignore this field if it is empty...
          else
              eval(['binned_site_info.' the_info_field_names{iInfo} '(i, :) = raster_site_info.' the_info_field_names{iInfo} ';']);   % might run into problems with this so above line could be more useful
          end
       
       end
       
   else
       binned_site_info = [];
   end
   

end


% save extra information about the bin_width, sampling_interval, etc.
binned_site_info.binning_parameters.raster_file_directory_name = raster_file_directory_name;
binned_site_info.binning_parameters.bin_width = bin_width;
binned_site_info.binning_parameters.sampling_interval = sampling_interval;

binned_site_info.binning_parameters.start_time = start_time;
binned_site_info.binning_parameters.end_time = end_time;

binned_site_info.binning_parameters.the_bin_start_times = the_bin_start_times;
binned_site_info.binning_parameters.the_bin_widths = the_bin_widths;


% if there is a field in the in the raster_site_info that has the time that was used to align the data
%  and if this time is the same for all raster files (as it should be) then save this time with the binning parameters
%  this time can then be used by plotting objects to align the data to specify what time is 0.
if isfield(binned_site_info, 'alignment_event_time')   
    if sum(abs(diff(binned_site_info.alignment_event_time))) == 0       
        binned_site_info.binning_parameters.alignment_event_time = raster_site_info.alignment_event_time(1);
    end
end



if (nargin < 5) && (length(bin_width) == 1) && (length(sampling_interval) == 1) 
    saved_binned_data_file_name = [save_prefix_name '_' num2str(bin_width) 'ms_bins_' num2str(sampling_interval) 'ms_sampled'];     
elseif (nargin > 4) && (length(bin_width) == 1) && (length(sampling_interval) == 1)  
    saved_binned_data_file_name = [save_prefix_name '_' num2str(bin_width) 'ms_bins_' num2str(sampling_interval) 'ms_sampled_' num2str(start_time) 'start_time_' num2str(end_time) 'end_time']; 
elseif (length(bin_width) == 1) && (length(sampling_interval) > 1)
    saved_binned_data_file_name = [save_prefix_name '_' num2str(bin_width) 'ms_bins_custom_sampling'];         
elseif (length(bin_width) > 1) && (length(sampling_interval) > 1)
    saved_binned_data_file_name = [save_prefix_name '_custom_bins_and_custom_sampling']; 
end


fprintf('\n')
disp(['  Saving the binned data to the file:  ' saved_binned_data_file_name])


%save(saved_binned_data_file_name, '-v7.3', 'binned_data', 'binned_labels', 'binned_site_info');
save(saved_binned_data_file_name, 'binned_data', 'binned_labels', 'binned_site_info');






function  binned_data = bin_one_site(raster_data, the_bin_start_times, the_bin_widths)  
% a helper function that bins the data for one site

  for c = 1:length(the_bin_start_times)      
      binned_data(:, c) = mean(raster_data(:, the_bin_start_times(c):(the_bin_start_times(c) + the_bin_widths(c) -1)), 2);            
  end

  
  

    
    

            
        
     
   


