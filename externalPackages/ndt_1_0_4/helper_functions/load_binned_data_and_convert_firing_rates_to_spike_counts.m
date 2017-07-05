function [binned_data_spike_counts binned_labels binned_site_info] = load_binned_data_and_convert_firing_rates_to_spike_counts(binned_data_file_name, bin_width)

% This function allows one to load binned-format data and convert the binned_data for each site into spike counts.  
%   The reason this function is needed is because paricular classifiers, such as the Poisson Naive Bayes classifier,
%   operate on integer data, however create_time_averaged_binned_data_from_raster data
%   creates firing rates that are not integers, thus by using this method, the firing rate data will be converted into integer spike counts).
%
%
% Inputs arguments:
%
%   1.  binned_data_file_name: the file name of data that is in binned format (that contains binned_data which has firing rate data)
%
% Optional input arguemnts: 
%   
%   bin_width: the bin width that each column of binned_data{iSite} was averaged over to create firing rates. If this argument is a vector, then the values
%       correspond to the bin widths of each column of binned_data{iSite}; if this value is a scalar then it is assumed that all bin widths are the same size.
%       If this argument is not given then this method assumes there is a variable binned_site_info.binning_parameters.bin_width that will be used as the bin width (this variable
%       will exist if the function create_binned_data_from_raster_data was used to create the binned data.
%       
% Output arguments: 
%
%   1. binned_data_spike_counts:  binned_data that has been converted into spike counts  
%   2. binned_labels:  the binned labels
%   3. binned_site_info: the binned site info
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



   load(binned_data_file_name) 


   if nargin < 2 
        bin_width = binned_site_info.binning_parameters.bin_width;   % will either be a number, or if custom bin widths are used, this will be a vector
   end


   c = 1;
   for iSite = 1:size(binned_data, 2)   % go through each site and convert it into 

       if isscalar(bin_width)
           binned_data_spike_counts{iSite} =  binned_data{iSite} .* bin_width;
       else
           binned_data_spike_counts{iSite} = binned_data{iSite} .* repmat(bin_width, size(binned_data{iSite}, 1), 1);
       end

       % check to make sure that binned_data_spike_counts only contains integers showing that the conversion worked correctly
       %if length(find(rem(binned_data_spike_counts{iSite}, 1) ~= 0))
       %if length(find(rem(binned_data_spike_counts{iSite}, 1) > 10^-12)) % changed it to this to deal with numerical precision issues  
       
       % changed to make this work with Octave (rounding errors can occur in both directions) 
       deviations_from_intergers_matrix = min(rem(binned_data_spike_counts{iSite}, 1), 1 - rem(binned_data_spike_counts{iSite}, 1));
       
       if length(find(deviations_from_intergers_matrix > 10^-12))
           sites_with_failed_conversions(c) = iSite;
           c = c + 1;
       else
           binned_data_spike_counts{iSite} = round(binned_data_spike_counts{iSite});  % round data to get rid of small numerical imprecisions
       end

   end


   if exist('sites_with_failed_conversions')
        warning(['The following sites did not correctly convert to spike counts: ' num2str(sites_with_failed_conversions)]);
   end

   
end
        
 