function [reduced_binned_data, reduced_binned_labels] = reduce_data_to_particular_trials(binned_data, binned_labels, labels_for_determing_which_trial_to_use, label_names_to_use)

% This function will reduce the binned_data and binned_labels so that only trials that match particular criteria are used.  
%  In particular, this function finds all trials in which labels_for_determing_which_trial_to_use match one of the labels given in 
%  label_names_to_use, and then returns reduced_binned_data and reduced_binned_labels that only have data and labels from these trials.
%
%  Input arguments:
%   
%   1. binned_data:  The data in binned-format that should be reduced to using only particular trials. 
%
%   2. binned_labels:  The labels in binned-format that should be reduced to only using particular trials.
%
%   3. labels_for_determing_which_trial_to_use:  One particular binned-format label (i.e., binned_labels.particular_labels)
%       that are used to determine which trials to use.  
%   
%   4. label_names_to_use:  A cell array of label numbers (or a vector of numbers if labels_for_determing_which_trial_to_use is a vector).  
%       All trials in the labels_for_determing_which_trial_to_use that contain any of the labels in label_names_to_use will be included in the reduced data.
%
%  Output arguments:
%
%   1. reduced_binned_data:  data in binned-format that includes only trials that match the specified criteria.  
%
%   2. reduced_binned_labels  labels in binned-format that includes only trials that match the specified criteria. 
%
%  
%  Example:
%
%   Suppose we have some data in binned-format in which the field binned_labels.stimulus_ID exists. Then running:
%
%   [reduced_binned_data, reduced_binned_labels] = reduce_data_to_particular_trials(binned_data, binned_labels, binned_labels.stimulus_ID, [1 3 5 7])
%
%   will return reduced_binned_data and reduced_binned_labels in which only the trials in which stimulus_ID has values of 1, 3, 5, or 7 are present
%   (and all other trials have been removed).  In particular, if binned_labels has other fields apart from stimulus_ID, then these other labels will also be reduced.


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



% if strings are given for the labels, convert them to numbers 
if iscell(labels_for_determing_which_trial_to_use{1})
    [labels_for_determing_which_trial_to_use string_to_number_mapping] = convert_label_strings_into_numbers(labels_for_determing_which_trial_to_use);
    label_names_to_use = find(ismember(string_to_number_mapping, label_names_to_use));
end




for iNeuron = 1:length(binned_data)
    
    
   % find trials in which the label matches particular criteria 
   the_inds = []; 
   for iIDsToUse = 1:length(label_names_to_use)   
       the_inds = union(the_inds, find(labels_for_determing_which_trial_to_use{iNeuron} == label_names_to_use(iIDsToUse)));
   end
       
   % 'NaN' is a label of a trial that should be used
   if sum(isnan(label_names_to_use)) > 0
       the_inds = union(the_inds, find(isnan(labels_for_determing_which_trial_to_use{iNeuron}) == 1));
   end
   
    
   reduced_binned_data{iNeuron} = binned_data{iNeuron}(the_inds, :);
   
   
   the_label_names = fieldnames(binned_labels);
   
   for iLabels = 1:length(the_label_names)
       
      eval(['reduced_binned_labels.' the_label_names{iLabels} '{iNeuron} = binned_labels.' the_label_names{iLabels} '{iNeuron}(the_inds);']);
            
   end
   
    
   
end


































