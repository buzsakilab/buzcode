classdef basic_DS < handle

%  basic_DS implements the basic functions of a
%  datasource (DS) object, namely, it takes binned data and labels, 
%  and through the get_data method, the object returns k 
%  leave-one-fold-out cross-validation splits of the data which can subsequently
%  be used to train and test a classifier.  The data in the population vectors is 
%  randomly selected from the larger binned data that is passed to the constructor
%  of this object.  This object can create both
%  pseduo-populations  (i.e., populations vector in which the recordings  were 
%  made on independent sessions but are treated as if they were recorded simultaneously)
%  and simultaneously populations in which neurons that were recorded together
%  always appear together in population vectors.  
%
% Like all DS objects, basic_DS implements the method get_data, which has the following form:  
%
%  [XTr_all_time_cv YTr_all XTe_all_time_cv YTe_all] = get_data(ds);   where:
%
%    a. XTr_all_time_cv{iTime}{iCV} = [num_features x num_training_points] is a
%         cell array that has the training data for all times and cross-validation splits
%    b. YTr_all = [num_training_point x 1] a vector of the training labels  (the same training labels are used at all times and CV splits)
%    c. XTe_all_time_cv{iTime}{iCV} = [num_features x num_test_points] is a
%         cell array that has the test data for all times and cross-validation splits;   
%    d. YTe_all = [num_test_point x 1] a vector has the test labels (the same test labels are used at all times and CV splits)
%
%
%  The constructor for this object has the form:
%
%   ds = basic_DS(binned_data_name, specific_binned_label_name, num_cv_splits, load_data_as_spike_counts), where:
%
%     a. binned_data_name: is string that has the name of a file that has data in binned-format, or is a cell array of binned-format binned_data      
%     b. specific_binned_labels_name: is a string containing a specific binned-format label name, or is a cell array/vector containing 
%           the specific binned names (i.e., binned_labels.specific_binned_label_name)
%     c. num_cv_splits = is a scalar indicating how many cross-validation splits there should be
%     d. load_data_as_spike_counts:  an optional flag that can be set that will cause the data to be converted to spike counts if set to an integer rather than 0
%           (the create_binned_data_from_raster_data function saves data as firing rates by default).  This flag is useful
%           when using the Poison Naive Bayes classifier that needs spike counts rather than firing rates. If this flag is not set, the default behavior
%           is to use firing rates.
%
%
%  The basic_DS also has the following properties that can be set:
%
%  1. create_simultaneously_recorded_populations (default = 0).  If the data from all sites
%      was recorded simultaneously, then setting this variable to 1 causes the 
%      function to return simultaneous populations rather than pseudo-populations 
%      (for this to work all sites in 'the_data' must have the trials in the same order).  
%      If this variable is set to 2, then the training set is pseudo-populations and the
%      test set is simultaneous populations.  This allows one to estimate I_diag, as
%      described by Averbeck, Latham and Pouget in 'Neural correlations, population coding
%      and computation', Nature Neuroscience, May, 2006.  I_diag is a measure that gives a 
%      sense of whether training on pseudo-populations leads to a the same decision rule as 
%      when training on simultaneous populations.  
%
%  2.  sample_sites_with_replacement (default = 0).  This variable specifies whether 
%        the sites should be sample with replacement - i.e., if the data is 
%        sampled with replacement, then some sites will be repeated within a single 
%        population vector.  This allows one to do a bootstrap estimate of variance
%        of the results if different sites from a larger population had been selected
%        while also ensuring that there is no overlapping data between the training
%        and test sets.  
%                    
%  3. num_times_to_repeat_each_label_per_cv_split (default = 1).  This variable 
%        specifies how many times each label should appear in each cross-validation split.
%        For example, if this value is set to k, this means that there will be k
%        population vectors from each class in each test set, and there will be 
%        k * (num_cv_splits - 1) population vectors for each class in each training set split.  
%     
%  4. label_names_to_use (default = [] meaning all unique label names in the_labels are used).  
%       This specifies which labels names (or numbers) to use, out of the unique label
%       names that are present in the the_labels cell array. If only a subset of labels are listed, 
%       then only population vectors that have the specified labels will be returned.  
%
%  5. num_resample_sites (default = -1, which means use all sites).  This variable specifies
%        how many sites should be randomly selected each time the get_data method is called.
%        For example, suppose length(the_data) = n, and num_resample_sites = k, then each
%        time get_data is called, k of the n sites would randomly be selected to be included
%        as features in the population vector.  
%      
%  6. sites_to_use (default = -1, which means select features from all sites).  This
%       variable allows one to only choose features from the sites listed in this vector
%       (i.e., features will only be randomly selected from the sites listed in this vector).
%    
%  7. sites_to_exclude (default = [], which means do not exclude any sites).  This allows          
%       one to not select features from particular sites (i.e., features will NOT be
%       selected from the sites listed in this vector).
%
%  8. time_periods_to_get_data_from (default = [], which means create one feature
%       for all times that are present in the_data{iSite} matrix).  This variable 
%       can be set to a cell array that contains vectors that specify which time bins 
%       to use as features from the_data. For examples, if time_periods_to_get_data_from = {[2 3], [4 5], [10 11]}
%       then there will be three time periods for XTr_all_time_cv and  XTe_all_time_cv 
%       (i.e., length(XTr_all_time_cv) = 3), and the population vectors for the 
%       time period will have 2 * num_resample_sites features, with the population 
%       vector for the first time period having data from each resample site from times
%       2 and 3 in the_data{iSite} matrix, etc..  
%
%  9. randomly_shuffle_labels_before_running (default = 0). If this variable is set to one
%       then the labels are randomly shuffled prior to the get_data method being called (thus all calls
%       to get_data return the same randomly shuffled labels). This method is useful for creating a 
%       null distribution to test whether a decoding result is above what one would expect by chance.  
%
%
%  This object also has two addition method which are:
%  
%  1.  the_properties = get_DS_properties(ds)
%       This method returns the main property values of the datasource.
%
%  2.  ds = set_specific_sites_to_use(ds, curr_resample_sites_to_use)
%        This method causes the get_data to use specific sites rather than
%        choosing sites randomly.  This method should really only be used by 
%        other datasources that are extending the functionality of basic_DS.
%
%
%  Note: this class is a subclass of the handle class, meaning that when this object is created a
%   reference to the object is returned.  Thus when fields of the object are changed a copy of the 
%   object does not need to be returned (by default matlab objects are passed by value). The
%   advantage of having this object inherit from the handle class is that if the object changes its
%   state within a method, a copy of the object does not need to be returned (this is particularly
%   useful for the randomly_permute_labels_before_running method so that the labels can be randomly
%   shuffled once prior to the get_data method being called, and the same shuffled labels will
%   be used throughout all subsequent calls to get_data, allowing one to create a full null distribution
%   by running the code multiple times).
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




  
properties 
    
    
    the_labels                    %  a cell array that contains vectors of labels that specify what occurred during each trial for all neurons the_data cell array  
    
    num_cv_splits                 %  how many cross-validation splits there should be 
    
    num_times_to_repeat_each_label_per_cv_split = 1;     %  how many of each unique label should be in each CV block
    
    label_names_to_use = [];     %  which set of labels names should be used (or which numbers should be used if the_labels{iSite} is a vector of numbers)

    num_resample_sites = -1;            %  how many sites should be used for each resample iteration - must be less than length(the_data)
    sample_sites_with_replacement = 0;  %  specify whether to sample neurons with replacement - if they are sampled with replacement, then some features will be repeated within a single data vector
                                        %    which reduces the number of neurons used (and thus usually lowers the decoding accuracy).  It should be noted that all data duplication appears in the same CV trials
                                        %    so there is no contamination with having repeatd data in different CV trials
    
    create_simultaneously_recorded_populations = 0;   % to use pseudo-populations or simultaneous populations (2 => that the training set is pseudo and test is simultaneous)
      
    sites_to_use = -1;                   %  a list of indices of which features (e.g., sites/neurons) to use in the the_data cell array 
    sites_to_exclude = [];               %  a list of features that should explicitly be excluded
    time_periods_to_get_data_from = [];  %  a cell array containing vectors that specify which time bins to use from the_data 
 
    
    % randomly shuffles the labels prior to the get_data method being called - which is useful for creating one point in a null distribution to check if decoding results are above what is expected by change.
    randomly_shuffle_labels_before_running = 0;    
       
    binned_site_info    %  a variable that can contain the binned_site_info (this will be automatically set the data in the constructor is loaded from a string that has a file name of data in binned-format.) 
                        %       this information is not used by the datasource, but it is returned by the get_properties method.   
   
    % excluding this option for now
    % set these if you want the datasource to return different random selection of labels (and data from those trials) each time
    % use_random_subset_of_k_labels_each_time_data_is_retrieved = -1;  % if this is set to k > 1, then a random subset of labels k labels will be chosen each resample iteration 
    %  10. use_random_subset_of_k_labels_each_time_data_is_retrieved (default = -1).  If this 
    %       variable is set to k > 1, then a random subset of labels k labels will be chosen each 
    %       time get_data is run (i.e., the different resample runs will use a different subset
    %       of k labels each time the get_data method is called).  

    
                                   
           
end



properties (GetAccess = 'public', SetAccess = 'private')
     
     the_data                      %  a cell array that contains all the data binned data in the format the_data{iNeuron}[num_trials x num_time_bins]      
        
     curr_resample_sites_to_use = -1;  % This specified which features should be used on a given random resample iteration.  
                                % This should really be randomly selected each time the get_data is called, but some rare cases
                                % it is useful to set to some specific features.  The method set_specific_sites_to_use that allows this variable
                                % to be set by outside calls.  
                                
      initialized = 0;
      label_names_to_label_numbers_mapping = [];
      
      
      data_loaded_as_spike_counts;   % records whether the data was loaded as spike counts (rather than firing rates)
      
end


methods

    
        %function ds = basic_DS
        %end
    
    
        % the constructor
        function ds = basic_DS(binned_data_name, specific_binned_label_name, num_cv_splits, load_data_as_spike_counts)
            
            if nargin < 4
                load_data_as_spike_counts = 0;
            end
            
            
            if load_data_as_spike_counts > 0
                
                if ~isstr(binned_data_name)
                      error('If the argument load_data_as_spike_counts is set to a value greater than 0, binned_data_name must be a string listing the name of a file that has data in binned format')  
                end
                
                [binned_data_spike_counts binned_labels binned_site_info] = load_binned_data_and_convert_firing_rates_to_spike_counts(binned_data_name);
                
                ds.the_data = binned_data_spike_counts;  
                ds.binned_site_info = binned_site_info;
            
                
            elseif isstr(binned_data_name)
                load(binned_data_name)
                ds.the_data = binned_data;
                ds.binned_site_info = binned_site_info;
            else
                ds.the_data = binned_data_name;
            end
                
            if isstr(specific_binned_label_name)
                ds.the_labels = eval(['binned_labels.' specific_binned_label_name]);
            else          
                ds.the_labels = specific_binned_label_name;
            end
                
            ds.num_cv_splits = num_cv_splits;     
            
            
            ds.data_loaded_as_spike_counts = load_data_as_spike_counts;   % might as well save this information too
            
                                                            
        end
        
        
 
        
        % This method allows one to set exact prespecified sites to get data from 
        %  rather than randomly selecting a set of sites from the larger population (it should rarely be used)
        function ds = set_specific_sites_to_use(ds, curr_resample_sites_to_use)
            ds.curr_resample_sites_to_use = curr_resample_sites_to_use;
        end
           
        
        
        %  This method returns the main property values of the datasource (could be useful for saving what parameters were used)
        function the_properties = get_DS_properties(ds)
            
            the_properties.num_cv_splits =  ds.num_cv_splits; 
            the_properties.num_times_to_repeat_each_label_per_cv_split = ds.num_times_to_repeat_each_label_per_cv_split;
            the_properties.sample_sites_with_replacement = ds.sample_sites_with_replacement;
            the_properties.num_resample_sites = ds.num_resample_sites;
            the_properties.create_simultaneously_recorded_populations = ds.create_simultaneously_recorded_populations;
            the_properties.sites_to_use = ds.sites_to_use;
            the_properties.sites_to_exclude = ds.sites_to_exclude;
            the_properties.time_periods_to_get_data_from = ds.time_periods_to_get_data_from;                         
            the_properties.randomly_shuffle_labels_before_running = ds.randomly_shuffle_labels_before_running;
            the_properties.binned_site_info = ds.binned_site_info;
            the_properties.data_loaded_as_spike_counts = ds.data_loaded_as_spike_counts;
            %the_properties.use_random_subset_of_k_labels_each_time_data_is_retrieved = ds.use_random_subset_of_k_labels_each_time_data_is_retrieved;

            
            % if haven't converted to ds.label_names_to_use to strings yet (b/c get_data has not yet been called), or if the_labels is numbers, just return input set by user
            if iscell(ds.label_names_to_use) || isempty(ds.label_names_to_label_numbers_mapping)
                 the_properties.label_names_to_use = ds.label_names_to_use;            
            else  % if code has already converted ds.label_names_to_use from strings to numbers convert them back to strings                                      
                for iName = 1:length(ds.label_names_to_use)
                    the_properties.label_names_to_use{iName} = ds.label_names_to_label_numbers_mapping{ds.label_names_to_use(iName)};  % is cell array so won't work 
                end
            end

        end
        
        
        
        function  [XTr_all_time_cv YTr_all XTe_all_time_cv YTe_all] = get_data(ds) 
        % The main DS function that returns training and test population vectors.  The outputs of this function are:
        %    
        %  1. XTr_all_time_cv{iTime}{iCV} = [num_features x num_training_points] is a
        %        cell array that has the training data for all times and cross-validation splits
        %
        %  2. YTr_all{iTime} = [num_training_point x 1] has the training labels
        %
        %  3. XTe_all_time_cv{iTime}{iCV} = [num_features x num_test_points] is a
        %       cell array that has the test data for all times and cross-validation splits   
        %
        %  4. YTe_all{iTime} = [num_test_point x 1] has the test labels
        %
 
        
            % initialize variables the first time ds.get_data is called
            if ds.initialized == 0
            
                disp('initializing basic_DS.get_data')
                
                                
                % if the_labels is a cell array of strings, convert the_labels into a vector of numbers
                if iscell(ds.the_labels{1})  % just checking the first site (assuming it will be the same for all other sites)
                    
                    ignore_case_of_strings = 0;  % for now, always respect the case of the strings used in the labels
                    
                    % % [specific_binned_labels_as_numbers string_to_number_mapping] = convert_label_strings_into_numbers(ds.the_labels, ignore_case_of_strings, ds.label_string_names_to_use);                    
                     
                    % doing it this way causes a consistent mapping from label strings to label numbers, and then the strings to be used are selected through setting ds.label_names_to_use
                    [specific_binned_labels_as_numbers string_to_number_mapping] = convert_label_strings_into_numbers(ds.the_labels, ignore_case_of_strings); 
                    
                    ds.the_labels = specific_binned_labels_as_numbers;
                    ds.label_names_to_label_numbers_mapping = string_to_number_mapping;
                    
                    
                    if ~isempty(ds.label_names_to_use)
                        
                        label_numbers_used = find(ismember(string_to_number_mapping, ds.label_names_to_use));
                        
                        % if a label_names_to_use name contains a string that is not one of the strings in the_labels, print an error message
                        inds_of_bad_string_to_use_names = find(~ismember(ds.label_names_to_use, string_to_number_mapping));
                        if  ~isempty(inds_of_bad_string_to_use_names)
                                                       
                            bad_string_names = '';                           
                            for iBadStringName = 1:length(inds_of_bad_string_to_use_names)
                                bad_string_names = [bad_string_names ' ' ds.label_names_to_use{inds_of_bad_string_to_use_names(iBadStringName)} ','];
                            end
                            
                            valid_string_names = '';
                            for iValidString = 1:length(string_to_number_mapping)
                                valid_string_names = [valid_string_names ' ' string_to_number_mapping{iValidString} ','];
                            end

                            error(['ds.label_string_names_to_use must be set to names in this list:' valid_string_names(1:end-1) '.  ' ...
                                'The following ds.label_string_names_to_use strings not in the list:' bad_string_names(1:end-1)]);                                 
                        end
                        
                        
                        ds.label_names_to_use =  label_numbers_used;  % convert given label names that should be used into numbers
                        
                    else
                        ds.label_names_to_use = 1:length(string_to_number_mapping);
                    end
                                        
                end
                
                
               
                
                % if ds.randomly_shuffle_labels_before_running == 1, randomly shuffle the labels the first time get_data is run
                if (ds.randomly_shuffle_labels_before_running == 1)
                    'randomly shuffling the labels'
                    
                    % added in NDT version 1.0.2 so that for simultaneously recorded populations the labels in all sites are shuffled the same way (previous version returned an error when shuffling labels on simultaneously recorded populations)
                    if ds.create_simultaneously_recorded_populations > 0   
                        
                        shuffled_labels = ds.the_labels{1}(randperm(length(ds.the_labels{1})));  % all sites should have the same labels, so will shuffle the labels for the first site only and will use this order for all sites                           
                        for iSite = 1:length(ds.the_labels)                     
                            ds.the_labels{iSite} =  shuffled_labels;
                        end
                     
                    % for non-simultaneously recorded datasets, shuffle each channel separately (same as NDT version 1.0.0)   
                    else
                        for iSite = 1:length(ds.the_labels)                 % will shuffle the labels from all sites, not just from those specified in sites_to_use     
                            ds.the_labels{iSite} =  ds.the_labels{iSite}(randperm(length(ds.the_labels{iSite})));
                        end
                    end
                    
                end
                
                
                
                % if using simultaneously recorded populations, convert data to a format that will make code run a little faster                                   
                if ds.create_simultaneously_recorded_populations == 0
                
                    if ~(iscell(ds.the_data))
                        the_data = ds.the_data;
                        the_labels = ds.the_labels;

                        for iSite = 1:size(the_data, 2)                       
                            curr_data{iSite} = squeeze(the_data(:, iSite, :));
                            curr_labels{iSite} = ds.the_labels;
                        end

                        ds.the_data = curr_data;
                        ds.the_labels = curr_labels;
                    end
                         
                        
                elseif ds.create_simultaneously_recorded_populations > 0
     
                    if iscell(ds.the_data)
                        simultaneous_labels_to_use = ds.the_labels{1};   % this should be ok, since all channels should have all the labels (regardless of whether a channel will ultimately be used)
                        the_data = ds.the_data;

                        for iSite = 1:length(the_data)                        

                                if sum(abs(simultaneous_labels_to_use - ds.the_labels{iSite})) ~= 0
                                    error('problem, all simultaneously recorded neurons should have the same labels')   
                                end

                                the_simultaneous_data(:, :, iSite) = the_data{iSite};                          
                        end

                        ds.the_data = permute(the_simultaneous_data, [1 3 2]);
                        ds.the_labels = simultaneous_labels_to_use;
                    end

                end
                
                
                if isempty(ds.time_periods_to_get_data_from)
                    
                    if iscell(ds.the_data)   % for pseudo-populations 
                        num_time_periods = size(ds.the_data{1}, 2);
                    else   % for simultaneous data 
                        num_time_periods = size(ds.the_data, 3);
                    end
                        
                    for i = 1:num_time_periods 
                        time_periods_to_use{i} = i;                
                    end
                    ds.time_periods_to_get_data_from = time_periods_to_use;         

                end
                
                
                
                % now that everything has been initialized, set inialized flag to 1     
                ds.initialized = 1;
              
                
            end   % end initialization
            
            
            
            
            % access to objects' fields in matlab is very slow (which is super pathetic), so to get the code to run faster I have use temporary copies of the data
            %  (hopefully matlab will fix this in the future)
            the_data = ds.the_data;                                                 
            the_labels = ds.the_labels;            
            num_cv_splits = ds.num_cv_splits;            
            num_times_to_repeat_each_label_per_cv_split =  ds.num_times_to_repeat_each_label_per_cv_split;                        
            curr_resample_sites_to_use = ds.curr_resample_sites_to_use;
            label_names_to_use = ds.label_names_to_use;  if size(label_names_to_use, 1) ~= 1, label_names_to_use = label_names_to_use'; end  % make sure labels numbers are in the correct orientation
            sites_to_use = ds.sites_to_use;
            sites_to_exclude = ds.sites_to_exclude;
            num_resample_sites = ds.num_resample_sites;
            sample_sites_with_replacement = ds.sample_sites_with_replacement;   
            create_simultaneously_recorded_populations = ds.create_simultaneously_recorded_populations;
            %use_random_subset_of_k_labels_each_time_data_is_retrieved = ds.use_random_subset_of_k_labels_each_time_data_is_retrieved;
                    
                    
            % a santy checks
            if isempty(sites_to_use) 
               error('sites_to_use can not be empty') 
            end

            
            % if sites_to_use is a number that is less than 0, use all sites
            if  (sites_to_use < 1)
                if create_simultaneously_recorded_populations > 0
                    sites_to_use = 1:size(the_data, 2);
                else
                    sites_to_use = 1:length(the_data);
                end    
            end
        
             
            if ~isempty(sites_to_exclude)      % can exclude specific neurons as well as specify which ones should be used
                sites_to_use = setdiff(sites_to_use, sites_to_exclude);
            end
            
            if isempty(label_names_to_use) ||  (isscalar(label_names_to_use) && label_names_to_use < 1)
               if create_simultaneously_recorded_populations > 0
                   label_names_to_use = unique(the_labels);
               else
                   label_names_to_use = unique(the_labels{sites_to_use(1)});   % use all the labels as the default value  (assuming that the first used neuron 1 has all the labels shown - which might not be a foolproof assumption)  % changed on 1/25/12               
               end                
            end
                         
            
            % more sanity checks
            if length(label_names_to_use) ~= length(unique(label_names_to_use))
               warning('some labels were listed twice in the field ds.label_names_to_use, (these duplicate enteries will be ignored)');
               label_names_to_use = unique(label_names_to_use);
            end
            
            
            % making sure create_simultaneously_recorded_populations is a valid argument            
            if (create_simultaneously_recorded_populations  > 2) || (create_simultaneously_recorded_populations  <  0)
                error('create_simultaneously_recorded_populations must be set to 0, 1 or 2');   
            end
            

            % if the number of resample neurons is not specified, use all neurons
            if num_resample_sites < 1  
                num_resample_sites = length(sites_to_use);
            end
            

            % code for randomly selecting k labels to use each time data is retrieved 
            %if use_random_subset_of_k_labels_each_time_data_is_retrieved > 1   % needs at least 2 labels for a classification problem to work
            %       rand_label = label_names_to_use(randperm(length(label_names_to_use)));
            %       label_names_to_use = rand_label(1:use_random_subset_of_k_labels_each_time_data_is_retrieved); 
            %end
            
            
            % if specific sites to be used have not been given (as should usually be the case), randomly select some sites to use from the larger population
            if  length(curr_resample_sites_to_use) == 1 && (curr_resample_sites_to_use < 1) %isempty(curr_resample_sites_to_use)

                if ~(sample_sites_with_replacement)   % only use each feature once in a population vector
                    curr_resample_sites_to_use = sites_to_use(randperm(length(sites_to_use)));
                    curr_resample_sites_to_use = sort(curr_resample_sites_to_use(1:num_resample_sites));  % sorting just for the heck of it

                else   % selecting random features with replacement (i.e., the same feature can be repeated multiple times in a population vector).    
                    initial_inds = ceil(rand(1, num_resample_sites) * num_resample_sites);  % can have multiple copies of the same feature within a population vector
                    curr_resample_sites_to_use = sort(sites_to_use(initial_inds));    
                end

            end

            
            % making code more robust in case transpose of curr_resample_sites_to_use is actually passed as an argument
            if (size(curr_resample_sites_to_use, 1) > 1)
                curr_resample_sites_to_use = curr_resample_sites_to_use';
            end
                  
           
           % pre-allocate memory
           all_data_point_labels = NaN .* ones(length(unique(label_names_to_use)) * num_cv_splits * num_times_to_repeat_each_label_per_cv_split, 1);           
           start_boostrap_ind = 1;  
            
           
           % make sure label_names_to_use is a row vector
           if size(label_names_to_use, 1) > 1
               label_names_to_use = label_names_to_use';
           end

         
           
           if  create_simultaneously_recorded_populations == 0

       
               % pre-allocate memory. 
                the_resample_data = NaN .* ones(length(unique(label_names_to_use)) * num_cv_splits * num_times_to_repeat_each_label_per_cv_split, length(curr_resample_sites_to_use), size(the_data{1}, 2)); 
          

                % if someone has changed the data or the labels after they have already set create_simultaneously_recorded_populations = 0, then the format of these variables needs to be converted back
                if ~(iscell(ds.the_data)) || ~(iscell(the_labels)) 
                    create_simultaneously_recorded_populations = 0;
                end
                 
                
               % create a 3 dimensional tensor the_resample_data that is [(num_labels * num_cv_slits * num_repeats_per_cv_label)   x num_neurons x num_time_bins] large               
               for iLabel = label_names_to_use 

                    cNeuron = 1;    
                    for iNeuron = unique(curr_resample_sites_to_use)  

                        % choose random trials to use for each label type
                        curr_trials_to_use = find(the_labels{iNeuron} == iLabel);  %find(ds.the_labels{iNeuron} == iLabel);
                        curr_trials_to_use = curr_trials_to_use(randperm(length(curr_trials_to_use)));

                        if length(curr_trials_to_use) < (num_cv_splits  * num_times_to_repeat_each_label_per_cv_split)    
                            error(['Requestion data from more trials of a given condition than has been recorded.  This is due to ' ...
                                '(ds.num_cv_splits * ds.num_times_to_repeat_each_label_per_cv_split) being greater than the number of times a given condition ' ...
                                'is present in the data (for at least one site). Make sure that only sites that have enough repetitions of each condition are used ' ...
                                '(this can be done by setting ds.site_to_use = find_sites_with_at_least_k_repeats_of_each_label(the_labels_to_use, num_cv_splits) )']);  
                            return;
                        else
                            curr_trials_to_use = curr_trials_to_use(1:(num_cv_splits  * num_times_to_repeat_each_label_per_cv_split));
                        end


                        % put everything into the correct number of CV splits
                        for iRepeats = 1:length(find(curr_resample_sites_to_use == iNeuron))                    
                            the_resample_data(start_boostrap_ind:(start_boostrap_ind + length(curr_trials_to_use) - 1), cNeuron, :) = the_data{iNeuron}(curr_trials_to_use, :);
                            cNeuron = cNeuron + 1;
                        end

                    end    % end iNeuron


                    all_data_point_labels(start_boostrap_ind:(start_boostrap_ind + length(curr_trials_to_use) - 1)) = iLabel .* ones(length(start_boostrap_ind:(start_boostrap_ind + length(curr_trials_to_use) - 1)), 1);

                    start_boostrap_ind = start_boostrap_ind + length(curr_trials_to_use);


                end   % end for iLabel  

                
                
            
            % if creating simultaneously recorded populations ...
           elseif  create_simultaneously_recorded_populations > 0
            
                               
               the_resample_data = NaN .* ones(length(unique(label_names_to_use)) * num_cv_splits * num_times_to_repeat_each_label_per_cv_split, length(curr_resample_sites_to_use), size(the_data, 3)); 
 
               the_data = the_data(:, curr_resample_sites_to_use, :);

               
                % choose (num_cv * num_repeats) random data points for each class            
                for iLabel = label_names_to_use

                    % choose random trials to use for each label type
                    curr_trials_to_use = find(the_labels == iLabel);
                    curr_trials_to_use = curr_trials_to_use(randperm(length(curr_trials_to_use)));

                    if length(curr_trials_to_use) < (num_cv_splits * num_times_to_repeat_each_label_per_cv_split)
                        error('problems:  asking for more trials of a given stimuli type then were recorded in the experiment');  % maybe this should be an error
                        return;
                    else
                        curr_trials_to_use = curr_trials_to_use(1:(num_cv_splits * num_times_to_repeat_each_label_per_cv_split));
                    end

                    the_resample_data(start_boostrap_ind:(start_boostrap_ind + length(curr_trials_to_use) - 1), :, :) = the_data(curr_trials_to_use, :, :);                
                    all_data_point_labels(start_boostrap_ind:(start_boostrap_ind + length(curr_trials_to_use) - 1)) = iLabel .* ones(length(start_boostrap_ind:(start_boostrap_ind + length(curr_trials_to_use) - 1)), 1);             
                    start_boostrap_ind = start_boostrap_ind + length(curr_trials_to_use);

                end
            
           end % end simultaneous populations
            

           
            clear the_data   % clear up some memory
            
            
            all_resample_data_inds = 1:size(the_resample_data, 1);
            
            time_periods_to_get_data_from = ds.time_periods_to_get_data_from;
            
            
            % convert the_resample_data into a training and splits 
            for iTimePeriod = 1:length(time_periods_to_get_data_from)
        
                curr_data = the_resample_data(:, :, time_periods_to_get_data_from{iTimePeriod});            
                the_resample_data_time = reshape(curr_data, [size(curr_data, 1) size(curr_data, 2) * size(curr_data, 3)]);

                
                if (create_simultaneously_recorded_populations > 1)  
                    the_site_ids = 1:size(curr_data, 2);
                    simul_to_pseudo_feature_to_siteID_mapping{iTimePeriod} = repmat(the_site_ids, [1 size(curr_data, 3)]);
                end
                
                
                cv_start_ind = 1;
                
                for iCV = 1:num_cv_splits
                   
                    curr_cv_inds = [];  %NaN .* ones(num_times_to_repeat_each_label_per_cv_split * length(ds.label_names_to_use), 1);
                    
                    for iNumRepeatsPerLabel = 1:num_times_to_repeat_each_label_per_cv_split
                        
                        curr_cv_inds  = [curr_cv_inds cv_start_ind:(num_cv_splits  * num_times_to_repeat_each_label_per_cv_split):size(the_resample_data, 1)];
                        cv_start_ind = cv_start_ind + 1;
                    end
                    
                    
                    
                    % these cells arrays contain the data for each CV splits separately,
                    % but don't need to create these, rather this function will just return the CV data divided into training and test sets
                    %    % % cross_validation_splits_all_time_periods{iCV}{iTimePeriod} = the_resample_data_time(curr_cv_inds, :);  % old get_resample_data7 format...
                    %    cross_validation_splits_all_time_periods{iTimePeriod}{iCV} = the_resample_data_time(curr_cv_inds, :)';  % should replace above with this soon
                    %    cross_validation_labels{iCV} = all_data_point_labels(curr_cv_inds);   % not sure why I need a separate one for each CV?

                   
                    % can actually just return these instead...
                    XTr_all_time_cv{iTimePeriod}{iCV} = the_resample_data_time(setdiff(all_resample_data_inds, curr_cv_inds), :)';  
                    XTe_all_time_cv{iTimePeriod}{iCV} = the_resample_data_time(curr_cv_inds, :)';
                    
                    %YTr_all{iTimePeriod} = all_data_point_labels(setdiff(all_resample_data_inds, curr_cv_inds));
                    %YTe_all{iTimePeriod} = all_data_point_labels(curr_cv_inds);
                    
                    % might as well return these as vectors (rather than cell arrays) since they are the same at all time periods
                    YTr_all = all_data_point_labels(setdiff(all_resample_data_inds, curr_cv_inds));  
                    YTe_all = all_data_point_labels(curr_cv_inds);
                    
                                 
                end                
                
            end
            
                       
            
            % If create_simultaneously_recorded_populations == 2, create pseudo-populations for training, and simultaneous data for testing.
            % This is useful for assessing I_diag as described by Averbeck, Latham and Pouget, Nature Neurosience, May 2006.
            if  create_simultaneously_recorded_populations == 2
                XTr_all_time_cv = turn_training_simultaneous_data_into_pseudo_populations(XTr_all_time_cv, YTr_all, simul_to_pseudo_feature_to_siteID_mapping);
            end
            
            % This is another type of sanity check on pseudo-populations, but really is no reason to use this, so I am not going to give this as an option for now
            % if  create_simultaneously_recorded_populations == 3
            %    [XTr_all_time_cv XTe_all_time_cv] = turn_all_simultaneous_data_into_pseudo_populations(XTr_all_time_cv, YTr_all, XTe_all_time_cv, YTe_all, simul_to_pseudo_feature_to_siteID_mapping)
            % end
  
            
        end   % end get_data
        
        
        
   end  % end methods
   
    
   
   
   
end % end class



