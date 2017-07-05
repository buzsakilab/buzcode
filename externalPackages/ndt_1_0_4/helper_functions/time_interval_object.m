classdef time_interval_object < handle
    
%  A helper object that creates time intervals based on bin widths, sampling intervals, etc..  These time intervals 
%    (vectors of times) can be used to plot decoding accuracies against or to get the latency when decoding results are above chance.  
%    The constructor of this object does not take any arguments (i.e., time_interval_obj = time_interval_object). 
%    The properties that dictate how the time intervals are constructed can either be set directly or 
%    can be set through structures using the set_parameters_with_binning_stucture and 
%    set_all_currently_unset_parameters_with_binning_structure methods.  Once these properties have been set
%    the method the_time_interval = time_interval_obj.get_time_interval can be called to return the time interval.  These
%    methods are described in more detail below.
%
%
%  The time_interval_object has the following properties that must be set (either directly or through the set parameter methods)
%       that will dictate how the time interval is constructed:
%   
%    1.  .sampling_interval:  specifies the sampling interval between successive data points.
%           If sampling_interval is a vector, then the values in this vector will specify the start times
%           for the bins (e.g., if sampling_interval = [1 200 500], then the first bin will start time 1,
%           the second bin will start at time 200, and the third bin will start at time 500).  
%
%    2.  .bin_width:  the size of each bin.  If bin_width is a vector that is the same size 
%           as sampling_interval, then different bin widths can be specified for different time periods 
%           (e.g., if bin_widths = [50 500, 100], then the first bin will be 50 ms in length, 
%           the second will be 500 ms, and the third over 100 ms).
%
%  The following properties can also be set to dictate how the interval is constructed (default property values are used if these are not specified):      
%   
%    3.  .start_time (default = 1). The time of the start of the time intervals.  
%        
%    4.  .end_time (default = 1000 * sampling_interval(1)).  The end time of the time interval.
%
%    5.  .alignment_event_time (default = 0).  The time that should be specified as time zero in the interval (e.g., if one wanted an interval
%           from -500 to X, then one could set start_time = 1, and the alignment_event_time = 501).  
% 
%    6.  .alignment_type (default = 2).  Specifies whether the time intervals should be aligned to the start of the time bins (alignment_type = 1), 
%           the middle of the time interval (alignment_type = 2), or the end of the time bins (alignment_type = 3).  
%
%
%  Apart from setting the properties of this object directly, one can set them via a structure using the two methods described below: 
%
%   1. time_interval_obj.set_parameters_with_binning_stucture(binning_parameters_structure)
%           The input argument to this method, binning_parameters_structure, should contain fields that have some of the same names as the 
%           6 properties listed above (e.g., .sampling_interval, .start_time, etc.). All fields in this structure that match the property names
%           in the time_interval object will be set.  In particular, if the create_binned_data_from_raster_data function is used to create binned-format data,
%           then a structure binned_site_info.binning_parameters will be created by this function which can be passed to this method to
%           recreate a time interval.  
%
%   2.  time_interval_obj.set_all_currently_unset_parameters_with_binning_structure(binning_parameters_structure)
%           This method is the similar to the time_interval_obj.set_parameters_with_binning_stucture method except that it only 
%           sets properties for this object that have not already been previously set (i.e., it does not overwrite values that already exist). 
%           This is useful if one has set a bunch of properties by hand and then one wants to set all the rest of the properties using a default structure
%           such as the binned_site_info.binning_parameters structure.
%
%
%  The main method to create the time interval is:
%
%    the_time_interval = time_interval_obj.get_time_interval .  The returned value the_time_interval is a vector of numbers
%       that specifies times when events occurred (e.g., it is useful for labeling the x-axis when plotting results, etc.).
%
%

    
%==========================================================================

%     This code is part of the Neural Decoding Toolbox.
%     Copyright (C) 2013 by Ethan Meyers (emeyers@mit.edu)
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
    
    
    properties 
              
                               
        alignment_type;  % 1=> align to the start of the bin, 2 => align to the middle of the bin, 3 => align to the end of the bin 
                         %  if this is not set the default used by get_time_interval is to align to the center of the bin (2)
        
        % basic properties of this object
        bin_width;
        sampling_interval;
        
        start_time;   % if not set, the default if not set is 1
        end_time;    % if not set, the default is to create an end time that is 1000 times sampling_interval(1)

        alignment_event_time;   % if not set, the default is 0

    end
    
    
    
    methods
        
        
        function time_interval_obj = time_interval_object()            
        end
        
        
        function set_parameters_with_binning_stucture(time_interval_obj, binning_parameters_structure)
            
            
            if isempty(binning_parameters_structure)
                return
            end
                
            
            %time_interval_obj_properties = fields(time_interval_obj);
            time_interval_obj_properties = fieldnames(time_interval_obj);  % changed this to fieldnames to make it compartibel with Octave

            
            % set all the properties in the time_interval_object to the properties specified by the 
            for iProperty = 1:numel(time_interval_obj_properties)
                
                if isfield(binning_parameters_structure, time_interval_obj_properties{iProperty})
                    eval(['time_interval_obj.' time_interval_obj_properties{iProperty} ' = binning_parameters_structure.'  time_interval_obj_properties{iProperty} ';'])
                end
                
            end

            
            % check if there are any parameters in binning_parameters_structure that are not valid fields in the time_interval_object, if so give a warning           
            %binning_parameters_structure_fields = fields(binning_parameters_structure);
             binning_parameters_structure_fields = fieldnames(binning_parameters_structure);  % changed this to fieldnames to make it compartibel with Octave

            
            fields_in_binning_structure_not_in_time_interval_obj = setdiff(binning_parameters_structure_fields, time_interval_obj_properties);
            
            % ignore field 'raster_file_directory_name', 'the_bin_start_times', 'the_bin_width'
            fields_in_binning_structure_not_in_time_interval_obj = setdiff(fields_in_binning_structure_not_in_time_interval_obj, {'raster_file_directory_name', 'the_bin_start_times', 'the_bin_widths'});   
            
            if ~isempty(fields_in_binning_structure_not_in_time_interval_obj)
               warning(['The following fields in binning_parameters_structure that are not being used by the time_interval_obj: ' ])
               disp(strvcat(fields_in_binning_structure_not_in_time_interval_obj))
            end
            
            
        end
        
        
        
        
        function set_all_currently_unset_parameters_with_binning_structure(time_interval_obj, binning_parameters_structure)
            
                                   
            if isempty(binning_parameters_structure)
                return
            end
            
            
            %time_interval_obj_properties = fields(time_interval_obj);
            time_interval_obj_properties = fieldnames(time_interval_obj);  % changed this to fieldnames to make it compatible with Octave

            
            % set all the properties in the time_interval_object to the properties specified by the 
            for iProperty = 1:numel(time_interval_obj_properties)
                
                current_property_is_not_yet_set = eval(['isempty(time_interval_obj.' time_interval_obj_properties{iProperty} ')']);
                
                if isfield(binning_parameters_structure, time_interval_obj_properties{iProperty}) && current_property_is_not_yet_set
                    eval(['time_interval_obj.' time_interval_obj_properties{iProperty} ' = binning_parameters_structure.'  time_interval_obj_properties{iProperty} ';'])
                end
                
            end

            
            % check if there are any parameters in binning_parameters_structure that are not valid fields in the time_interval_object, if so give a warning           
            %binning_parameters_structure_fields = fields(binning_parameters_structure);
            binning_parameters_structure_fields = fieldnames(binning_parameters_structure);  % changed this to fieldnames to make it compatible with Octave
            
            fields_in_binning_structure_not_in_time_interval_obj = setdiff(binning_parameters_structure_fields, time_interval_obj_properties);

            % ignore field 'raster_file_directory_name', 'the_bin_start_times', 'the_bin_width'
            fields_in_binning_structure_not_in_time_interval_obj = setdiff(fields_in_binning_structure_not_in_time_interval_obj, {'raster_file_directory_name', 'the_bin_start_times', 'the_bin_widths'});   
            
            if ~isempty(fields_in_binning_structure_not_in_time_interval_obj)
               warning(['The following fields in binning_parameters_structure that are not being used by the time_interval_obj: ' ])
               disp(strvcat(fields_in_binning_structure_not_in_time_interval_obj))
            end
            
        end
        
        
        
        function the_time_interval = get_time_interval(time_interval_obj)            
        %  Use bin_width, sampling_interval, start_time, end_time to create a time interval 
        %  (ignoring binning_parameters.the_bin_start_times binning_parameters.the_bin_widths that are created by create_binned_data_from_raster_data
        %   b/c this create an ambiguity between which parameters to use to create the time interval, and because these are harder to set manually then
        %   setting the bin_width, sampling_interval, etc.).  
        
        
             if isempty(time_interval_obj.start_time)
                   start_time = 1;
             else
                   start_time = time_interval_obj.start_time;
             end


             if isempty(time_interval_obj.end_time)
                   %warning('Arbitarily creating a time interval that has 1000 time points.  One should set the time_interval_obj.end_time property to control how many time points are created.')
                   end_time = time_interval_obj.sampling_interval(1) .* 1000;   
             else
                   end_time = time_interval_obj.end_time;
             end


            % if a single bin width and step size have been specified, then create binned data that averaged data over bin_width sized bins, sampled at sampling_interval intervals
            if (length(time_interval_obj.bin_width) == 1) && (length(time_interval_obj.sampling_interval) == 1)  

                the_bin_start_times = start_time:time_interval_obj.sampling_interval:(end_time - time_interval_obj.bin_width  + 1);
                the_bin_widths = time_interval_obj.bin_width .* ones(size(the_bin_start_times));  


            % if multiple step size have been specified indicating bin start times, but only a single bin width is given, 
            %   then create binned_interval using bin_width sized bins, sampled using the start bin times given by sampling_interval
            elseif (length(time_interval_obj.bin_width) == 1) && (length(time_interval_obj.sampling_interval) > 1) 

                the_bin_start_times = time_interval_obj.sampling_interval;
                the_bin_widths = time_interval_obj.bin_width .* ones(size(the_bin_start_times)); 

            % if custom interval have been specified, use them
            elseif (length(time_interval_obj.bin_width) > 1) && (length(time_interval_obj.sampling_interval) > 1)     

                if length(time_interval_obj.bin_width) ~= length(time_interval_obj.sampling_interval), error('If a series of bin widths are specified by having bin_width be a vector, and a series of start binning times are specified by having sampling_interval be a vector, then the length(bin_width) must equal length(sampling_interval)'); end
                the_bin_start_times = time_interval_obj.sampling_interval;
                the_bin_widths = time_interval_obj.bin_width;

            end



            % sanity check that the_bin_start_times and the_bin_widths
            if numel(the_bin_start_times) ~= numel(the_bin_widths)
                error('the_bin_start_times and the_bin_widths must be the same size in order to create a time interval')
            end



            % get alignment type (times should be aligned to the end, center or beginning of the bins                 
            alignment_type = time_interval_obj.alignment_type;
            if isempty(alignment_type)    % align to the center of the bins by default
                alignment_type = 2;    
            end


            if alignment_type == 1
                the_time_interval = the_bin_start_times;   % align to the beginning of the bins
            elseif alignment_type == 2
                the_time_interval = the_bin_start_times + the_bin_widths./2;  % align to the center of the bins
            elseif alignment_type == 3
                the_time_interval = the_bin_start_times + the_bin_widths;  % align to the end of the bins
            else
                error('alignment_type must be set to 1, 2, or 3 corresponding to aligning the time interval to the end, middle or beginning of the time bins')
            end


            % get the time when the data should be aligned to
            alignment_event_time = time_interval_obj.alignment_event_time;
            if isempty(alignment_event_time)   % have 0 be the alignment event time by default
                alignment_event_time = 0;
            end 

            the_time_interval = the_time_interval - alignment_event_time;
        
        end
        

        
    end   % end methods
    
    
    
end   % end class