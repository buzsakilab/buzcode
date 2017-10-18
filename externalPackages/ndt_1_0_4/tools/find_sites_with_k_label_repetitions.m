function [inds_of_sites_with_at_least_k_repeats min_num_repeats_all_sites num_repeats_matrix label_names_used] = find_sites_with_k_label_repetitions(the_labels, k, label_names_to_use)

% This function takes in labels in binned label format, and an integer k, and returns the indices for all sites (e.g. neurons) 
%   that have at least k presentations of each condition.  The arguments to this function are:
%
%   1. the_labels: specific labels in binned-format that should be used (e.g., binned_labels.the_labels_to_use).  
%      
%   2. k: an integer specifying that each site returned should have at least k repetitions of each condition.
%
%  Optional input arguments:
%
%   3. label_names_to_use:  This specifies what label names (or numbers) to use (e.g., if the_labels contains labels the strings in the set {'red', 'green', 'blue'} 
%       but one only wants to know which sites have k repeats of 'red' and 'green' trials, then setting this to {'red', 'green'} will accomplish this goal). 
%       If this argument is not specified, then any label that was presented to any site will be used.
%
%
%  Outputs: 
%
%   1. inds_of_sites_with_at_least_k_repeats:  The indices of sites that have at least k repetitions of each condition.
%
%   2. min_num_repeats_all_sites:  This parameter lists for all sites the number of repetitions present for the label that
%       has the minium number of repetitions.
%
%   3. num_repeats_matrix:  A [num_sites x num_labels] matrix that specifies for each the number of repetitions of each condition
%       for each site (this variable could be useful for determining if particular conditions should be excluded based on whether
%       they were presented only a few times to many sites).
%
%   4. label_names_used: the names of the of the lables that were used (this is equivalent to label_names_to_use if this was passed
%       as an input argument). 
%
%  Example:  
%
%   Suppose we had an experiment in which a number of different stimuli were shown when recordings were made from a number
%       of different sites, and this information was contained in the variable binned_labels.stimulus_ID. The following command would
%       find all sites in which each stimulus condition was presented at least 20 times:
%
%       inds_of_sites_with_at_least_k_repeats = find_sites_with_at_least_k_repeats_of_each_label(binned_labels.stimulus_ID, 20)
%
%   When one is first starting to analyze a new dataset, one can also use this function to assess how many times each condition has
%       been presented to each site in order to determine how many CV splits to use.  Examining the variable min_num_repeats_all_sites
%       could be useful for this purpose, or one could run the following command:
%
%       for k = 0:60
%           inds_of_sites_with_at_least_k_repeats = find_sites_with_at_least_k_repeats_of_each_label(binned_labels.stimulus_ID, k);
%           num_sites_with_k_repeats(k + 1) = length(inds_of_sites_with_at_least_k_repeats);
%       end
%
%       The variable num_sites_with_k_repeats(i) indicates how many sites have at least i - 1 repetitions, i.e., num_sites_with_k_repeats(1) 
%       gives the total number of sites, num_sites_with_k_repeats(2) gives how many sites have at least one presentation of each stimulus, etc.
%       (note that 2 repetitions is the minimum needed to do a decoding analyses, although to get reasonable results usually more than 2 repetitions
%       are needed).



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




ignore_case_of_strings = 0;  % could make this an input argument, but not going to bother for now


if nargin < 3
    label_names_to_use = [];
    %label_names_to_use = {};  % switched this to make the code compatible with Octavem but problematic because now won't work with numbers as stimulus labels...
    for iSite = 1:length(the_labels)
        
        if ignore_case_of_strings == 1
            the_labels{iSite} = lower(the_labels{iSite});
        end
        
        label_names_to_use = union(label_names_to_use, the_labels{iSite});

    end
end



% create a [num_sites x num_labels_to_use] matrix called num_repeats_matrix that has how many times each each label was repeated for each site
for iSite = 1:length(the_labels)                
    for iLabel = 1:length(label_names_to_use)         
        if iscell(label_names_to_use)
            num_repeats_matrix(iSite, iLabel) = sum(ismember(the_labels{iSite}, label_names_to_use{iLabel}));
        else
            num_repeats_matrix(iSite, iLabel) = length(find(the_labels{iSite} == label_names_to_use(iLabel)));
        end
    end
end



if length(label_names_to_use) == 1
    min_num_repeats_all_sites = num_repeats_matrix';
else
    min_num_repeats_all_sites = min(num_repeats_matrix');
end



inds_of_sites_with_at_least_k_repeats = find(min_num_repeats_all_sites >= k);


% return the name of the labels that were used (useful if label_names_to_use was not an input parameter)
label_names_used = label_names_to_use;








