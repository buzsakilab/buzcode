function [specific_binned_labels_as_numbers string_to_number_mapping] = convert_label_strings_into_numbers(specific_binned_labels_as_strings, ignore_case_of_strings, string_to_number_mapping)

% This function takes specific binned labels (i.e., binned_labels.specific_binned_labels) where labels for each site are cell
%   arrays of strings and converts them into specific binned labels that are vectors of numbers (i.e., each unique string has been 
%   replaced by a unique number) These vector of numbers can then be used by other components of the toolbox (such as classifiers) to do
%   decoding analyses.
%
%  Input arguments:
%  
%    1. specific_binned_labels_as_strings: this is binned format labels that consist of strings (binned_labels.specific_labels), 
%        i.e., this is a cell array containing cell arrays of string for each site, where binned_labels.specific_labels{iSite} = a cell array of strings.
%        The cell arrays for each site will be converted into vectors of numbers, i.e., specific_binned_labels_as_numbers{iSite} = a vector of numbers.
%
%
%  Optional input arguments: 
%
%    2. ignore_case_of_strings (default = 0): if this is set to 1, the case of strings in 
%         specific_binned_labels_as_strings{iSite} will be ignored (e.g., 'LABELName' and 'Labelname' will
%         both be assigned to the same number in specific_binned_labels_as_numbers{iSite})
%
%    3. string_to_number_mapping: a cell array of strings where the order of the names in the cell array indicates
%         which number a given string will be assigned to (i.e., the string in string_to_number_mapping{3} will be given
%         be given the number 3 in specific_binned_labels_as_numbers{iSite}). If not all the unique names in specific_binned_labels_as_strings
%         are included in the string_to_number_mapping, then any missing name will be mapped to NaN.  If this argument is not given (or is an emptry matrix)
%         then strings will be assigned to numbers based on their alphabetic order (where capitalized letter occur before lower case letters if ignore_case = 0).  
%
%
%  Output arguments
%
%    1. specific_binned_labels_as_numbers: a cell array for each site, in which all the labels for each class are given in a vector of numbers 
%        (i.e., the cell array of strings for each site, specific_binned_labels_as_strings{iSite}, has now been converted into a vector of numbers).  
%
%    2. string_to_number_mapping: the order in which strings have been mapped into numbers, i.e., the string string_to_number_mapping{7} will have been
%        mapped to the number 7 in the specific_binned_labels_as_numbers{iSite} vector.
%


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




% default is to give strings will different capitalization unique numbers
if nargin < 2 || isempty(ignore_case_of_strings)
    ignore_case_of_strings = 0;
end


% get all the unique label names
all_unique_names = [];
for iSite = 1:length(specific_binned_labels_as_strings)
    
    if ignore_case_of_strings == 1
        specific_binned_labels_as_strings{iSite} = lower(specific_binned_labels_as_strings{iSite});
    end
        
    curr_unique_names = unique(specific_binned_labels_as_strings{iSite});    
    all_unique_names = union(all_unique_names, curr_unique_names);
       
end



% if string_to_number_mapping has been given as an input argument then use that mapping, otherwise use the mapping creating by finding all the unique labels
if ~exist('string_to_number_mapping') || isempty(string_to_number_mapping)
    string_to_number_mapping = all_unique_names;
else
     
    
   if length(unique(string_to_number_mapping)) ~= length(string_to_number_mapping)
      error('can not give input argument string_to_number_mapping that has redundant label names in the list') 
   end
    
    if ignore_case_of_strings == 1
        string_to_number_mapping = lower(string_to_number_mapping);
    end
    
    if length(string_to_number_mapping) ~= length(unique(string_to_number_mapping))
       error('can not have ignore_case_of_strings == 1 and have some of the same upper and lower case label names in the input argument string_to_number_mapping')
    end
        
end




% convert the specific binned labels as strings into specific binned labels as numbers
for iSite = 1:length(specific_binned_labels_as_strings)
       
    specific_binned_labels_as_numbers{iSite} = NaN .* ones(length(specific_binned_labels_as_strings{iSite}), 1); 
    
    for iLabel = 1:length(string_to_number_mapping)        
        specific_binned_labels_as_numbers{iSite}(ismember(specific_binned_labels_as_strings{iSite}, string_to_number_mapping{iLabel})) = iLabel;
    end
    
end



% if string_to_number_mapping has been given as an input (i.e., nargin = 3), 
%  and there are some names in string_to_number_mapping that are not in the_labels, give an error message
if length(setdiff(string_to_number_mapping, all_unique_names)) > 0
    
     % if a label_names_to_use name contains a string that is not one of the strings in the_labels, print an error message
    inds_of_bad_string_to_use_names = find(~ismember(string_to_number_mapping, all_unique_names));
    if  ~isempty(inds_of_bad_string_to_use_names)

        bad_string_names = '';                           
        for iBadStringName = 1:length(inds_of_bad_string_to_use_names)
            bad_string_names = [bad_string_names ' ' string_to_number_mapping{inds_of_bad_string_to_use_names(iBadStringName)} ','];
        end

        valid_string_names = '';
        for iValidString = 1:length(string_to_number_mapping)
            valid_string_names = [valid_string_names ' ' all_unique_names{iValidString} ','];
        end

        error(['All possible label names in binned format are:' valid_string_names(1:end-1) '.  ' ...
            'The following labels given are not in the list:' bad_string_names(1:end-1)]);                                 
    end
    
    
end









