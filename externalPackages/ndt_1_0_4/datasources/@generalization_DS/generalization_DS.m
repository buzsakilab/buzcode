classdef generalization_DS < handle
    
%  This datasource object (DS) allows one to train a classifier on a specific 
%  set of labels, and then test the classifier on a different set of 
%  labels - which enables one to evaluate how similar neural representations 
%  are across different but related conditions (i.e., does training on one set of
%  conditions generalization to a different but related set of conditions?). This datasource
%  is a subclass of the handle class (i.e., it has a persistent state) and 
%  contains a basic_DS where it gets most of its functionality.  
%
%  The constructor for this datasource has the same arguments
%  as basic_DS, plus two additional arguments 'the_training_label_names' 
%  and 'the_test_label_names' i.e., the constructor has the form:  
%
%     ds = generalization_DS(binned_data_name, specific_binned_label_name, num_cv_splits, the_training_label_names, the_test_label_names, load_data_as_spike_counts)
%   
%      the_training_label_names and the_test_label_names are cell arrays that
%        specify which labels should belong to which class, with the first element
%        of these cells arrays specifying the training/test labels for the first class
%        the second element of the cell array specifies which labels belong to 
%        the second class, etc..  For example, suppose one was interested in testing
%        position invariance, and had done an experiment in which data was recorded 
%        while 7 different objects were shown at three different locations.  If the
%        labels for the 7 objects at the first location had labels 'obj1_loc1', 'obj2_loc1', ..., 'obj7_loc1',
%        at the second location were 'obj1_loc2', 'obj2_loc2', ..., 'obj7_loc2',
%        and at the third location were 'obj1_loc3', 'obj2_loc3', ..., 'obj7_loc3',
%        then one could do a test of position invariance by setting the_training_label_names{1} = {'obj1_loc1},
%        setting the_training_label_names{2} = {'obj2_loc1'}, ...,  the_training_label_names{7} = {'obj7_loc1},
%        and setting the the_test_label_names{1} = {'obj1_loc2', 'obj1_loc3'}, 
%        the_test_label_names{2} = {'obj2_loc2', 'obj2_loc3'}, ..., the_test_label_names{7} = {'obj7_loc2', 'obj7_loc3'}. 
%        The object is able to test such generalization from training on one set of labels and
%        testing on a different set of labels by remapping the training label numbers to the
%        index number in the_training_label_names cell array, and remapping the 
%        test label numbers with the the index number into the the_test_label_names cell array.  
%
%  This object has all of the same properites as the basic_DS object (except that label_names_to_use which 
%    has been replaced by the the_training_label_names, and the_test_label_names properties).  There
%    is also an additional property that can be set for this object which is:  
%  
%       use_unique_data_in_each_CV_split  (default value is 0).
%  
%       When this argument is set to 0, the get_data method returns the normal leave one split
%        out training and test data sets (i.e., the training set consists of 
%        (num_cv_splits - 1) splits of the data and the test set consists of 1 split of the data).
%        The data in the training still comes from different splits as the data in the 
%        test set, thus one can have some of the same labels in the both 
%        the_training_label_names and in the_test_label_names (in fact, if ones sets
%        the_test_label_names = the_training_label_names, then the get_data
%        method will be the same as the basic_DS get_data method).  However, 
%        if use_unique_data_in_each_CV_split = 1, then each training 
%        and test set will consist data from only split, and thus each cross-validation
%        run is essentially like running an independent decoding experiment.  
%        In this case the_training_label_names and the_test_label_names must not contain 
%        any of the same labels (otherwise, they would be copies of the same data 
%        which would violate the fact that the training and the test set must not have any of the same data).  


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
    
    
    the_basic_DS = [];      % the basic_DS that will be used to give this object most of its functionality
    
    the_training_label_names = [];   % a cell array specifying which label names (or numbers) should be used for training the classifier
                                       %   i.e., the_training_label_names{1} = {class_1_training_names}; the_training_label_names{2} = {class_2_training_names}; etc.
    the_test_label_names = [];       % a cell array specifying which label names (or numbers) should be used for testing the classifier
                                       %   i.e., the_test_label_names{1} = {class_1_test_names}; the_test_label_names{2} = {class_2_test_names}; etc.
                                       
    use_unique_data_in_each_CV_split = 0;  % if this is set to 1 then each CV splits has unique data, i.e.,
                                                         %  each CV split has the amount of data that is typically only in the test set, 
                                                         %  and every CV training set does not consist of (num_cv_splits -1) * num_labels data points
                                                         %  but instead consists of length(cell2mat(the_test_label_names)) training points.
                                                         
                                                         
    % some properties of basic_DS that will also be available in generalization_DS by setting the basic_DS properties
                                                         
    num_times_to_repeat_each_label_per_cv_split = 1;     %  how many of each unique label should be in each CV block
    
    sample_sites_with_replacement = 0;  %  specify whether to sample neurons with replacement - if they are sampled with replacement, then some features will be repeated within a single point data vector
    num_resample_sites = -1;            %  how many sites should be used for each resample iteration - must be less than length(the_data)
    
    create_simultaneously_recorded_populations = 0;   % to use pseudo-populations or simultaneous populations (2 => that the training set is pseudo and test is simultaneous)
      
    sites_to_use = -1;                   %  a list of indices of which sites (neurons) to use in the the_data cell array 
    sites_to_exclude = [];               %  a list of features that should explicitly be excluded
    time_periods_to_get_data_from = [];  %  a cell array containing vectors that specify which time bins to use from the_data 

                                        
    % randomly shuffles the labels prior to the get_data method being called - which is useful for creating one point in a null distribution to check if decoding results are above what is expected by change.
    randomly_shuffle_labels_before_running = 0;                                                        
                                                         
                                  
end



properties (GetAccess = 'public', SetAccess = 'private')
 
      initialized = 0;
            
      the_training_label_numbers = [];   % numbers that the_training_label_names were mapping on to
      the_test_label_numbers = [];       % numbers that the_test_label_names were mapping on to
end

    

methods 

    
    function ds = generalization_DS(binned_data_name, specific_binned_label_name, num_cv_splits, the_training_label_names, the_test_label_names, load_data_as_spike_counts)
 
       if nargin < 6
           load_data_as_spike_counts = 0;
       end
        
       
       ds.the_basic_DS = basic_DS(binned_data_name, specific_binned_label_name, num_cv_splits, load_data_as_spike_counts);   % set properties using parent constructor
                      
       ds.the_training_label_names = the_training_label_names;   % set properties that are unique to this object
       ds.the_test_label_names = the_test_label_names;
          
       % get all the unique labels
       if length(the_training_label_names) ~= length(the_test_label_names)
            error('The cell arrays the_training_label_names and the_test_label_names must have the same number of cells (with the names/numbers in cell i containing which labels belong to class i)');
       end
       

    end
    

       
    function the_properties = get_DS_properties(ds)    
    
        the_properties = ds.the_basic_DS.get_DS_properties;
        
        the_properties.the_training_label_names = ds.the_training_label_names;
        the_properties.the_test_label_names  = ds.the_test_label_names; 
        
        the_properties.use_unique_data_in_each_CV_split = ds.use_unique_data_in_each_CV_split;
        
    end
    
    
    
    function  [XTr_all_time_cv YTr_all XTe_all_time_cv YTe_all] = get_data(ds)
    %function  [XTr_all_time_cv YTr_all XTe_all_time_cv YTe_all ADDITIONAL_DATASOURCE_INFO] = get_data(ds)
        

        
        % initialize things here....
        if ds.initialized == 0   
            
            
            % creating separate variables for this since calling a field of an object in Matlab is slow 
            the_training_label_names = ds.the_training_label_names;
            the_test_label_names = ds.the_test_label_names;
        
 
            ds.the_basic_DS.num_times_to_repeat_each_label_per_cv_split = ds.num_times_to_repeat_each_label_per_cv_split;

            ds.the_basic_DS.sample_sites_with_replacement = ds.sample_sites_with_replacement ;
            ds.the_basic_DS.num_resample_sites = ds.num_resample_sites;
            
            ds.the_basic_DS.create_simultaneously_recorded_populations = ds.create_simultaneously_recorded_populations;

            ds.the_basic_DS.sites_to_use = ds.sites_to_use;
            ds.the_basic_DS.sites_to_exclude = ds.sites_to_exclude;
            ds.the_basic_DS.time_periods_to_get_data_from = ds.time_periods_to_get_data_from;

            ds.the_basic_DS.randomly_shuffle_labels_before_running = ds.randomly_shuffle_labels_before_running;                                                  

    
            
            % if the_labels are strings, convert these names into label numbers...
            if iscell(ds.the_basic_DS.the_labels{1})
                
                
                % put all the training and test names into one cell array called all_training_and_test_names
                cTrainingNames = 0;
                cTestNames = 0;
                for iClass = 1:length(the_training_label_names)
                    for iTrain = 1:size(the_training_label_names{iClass}, 2)
                        cTrainingNames = cTrainingNames + 1;
                        all_training_names{cTrainingNames} = the_training_label_names{iClass}{iTrain};
                    end
                    
                    for iTest = 1:size(the_test_label_names{iClass}, 2)
                        cTestNames = cTestNames + 1;
                        all_test_names{cTestNames} = the_test_label_names{iClass}{iTest}; 
                    end                  
                end
                    
                
               all_training_and_test_names = union(all_training_names, all_test_names);

               
               % sanity checka that the same label is not in multiple training classes (and the the same label is not in multiple test classes)
               if length(unique(all_training_names)) ~= length(all_training_names)
                    warning('some of the same strings in the_training_label_names are in multiple classes');
               end  
               if length(unique(all_test_names)) ~= length(all_test_names)  % same sanity check for the_test_label_names
                    warning('some of the same strings in the_test_label_names are in multiple classes');
               end
               
               

               % convert the_labels into numebrs
               ignore_case_of_strings = 0;  % for now, always respect the case of the strings used in the labels

                
                % by passing all_training_and_test_names as the 3rd argument convert_label_strings_into_numbers this keeps the mapping from the_test_label_names to the_test_label_numbers correct                
                % an erorr in convert_label_strings_into_numbers will be thrown if some of the training or test label names are strings that are not in binned labels the_labels
                the_labels_as_numbers = convert_label_strings_into_numbers(ds.the_basic_DS.the_labels, ignore_case_of_strings, all_training_and_test_names); 

                ds.the_basic_DS.the_labels = the_labels_as_numbers;
                 
                
                % remap the_training_label_names and the_training_label_names into numbers
                cTrainingNames = 0;
                cTestNames = 0;
                for iClass = 1:length(the_training_label_names)
                    for iTrain = 1:size(the_training_label_names{iClass}, 2)
                        cTrainingNames = cTrainingNames + 1;
                        the_training_label_numbers{iClass}(iTrain) = find(ismember(all_training_and_test_names, the_training_label_names{iClass}{iTrain}));
                    end

                    for iTest = 1:size(the_test_label_names{iClass}, 2)
                        cTestNames = cTestNames + 1;
                        the_test_label_numbers{iClass}(iTest) = find(ismember(all_training_and_test_names, the_test_label_names{iClass}{iTest}));
                    end      
                end

                
               %  sanity check to make sure that one is not training with label l in class k and then has test label l in class j
               some_of_the_same_labels_are_in_multiple_classes = 0;
               for iClass = 1:length(the_training_label_names)                   
                    temp_test_label_numbers = the_test_label_numbers;
                    temp_test_label_numbers{iClass} = [];
                    temp_test_label_numbers = cell2mat(temp_test_label_numbers);
                    if ~isempty(intersect(temp_test_label_numbers, the_test_label_numbers{iClass}))
                        some_of_the_same_labels_are_in_multiple_classes = 1;
                    end                    
               end  
               
               if  some_of_the_same_labels_are_in_multiple_classes == 1
                    warning('some labels that are in training class k, are in a different test class j')
               end
               
                      
               ds.the_training_label_numbers = the_training_label_numbers;
               ds.the_test_label_numbers = the_test_label_numbers;
                 
                              
            else  % if numbers for labels have been specified instead of names, just use those numbers
                
               ds.the_training_label_numbers = ds.the_training_label_names;
               ds.the_test_label_numbers = ds.the_test_label_names;
               
            end
            

            % only use labels that are listed in the_training_label_numbers and the_test_label_numbers
            the_training_nums = cell2mat(ds.the_training_label_numbers);
            the_test_nums = cell2mat(ds.the_test_label_numbers);
            label_numbers_to_use = unique([the_training_nums(:); the_test_nums(:)]);
            ds.the_basic_DS.label_names_to_use = label_numbers_to_use;
 
            
            ds.initialized = 1;         
        end
        
        
                    
         % creating separate variables for this since calling a field of an object in Matlab is slow 
         the_training_label_numbers = ds.the_training_label_numbers;
         the_test_label_numbers = ds.the_test_label_numbers;

                
        
        %[XTr_all_time_cv YTr_all XTe_all_time_cv YTe_all] = get_data@basic_DS(ds);  % old version where this object was a subclass of basic_DS
         [XTr_all_time_cv YTr_all XTe_all_time_cv YTe_all] = ds.the_basic_DS.get_data;  % this might be much slower than inheriting from basic_DS :(
        
         
        if ds.use_unique_data_in_each_CV_split == 1
            
            
            % if running each training and test set separately, there can not be any overlap between the labels listed in the training and test sets 
            %(otherwise there could be some of the same data in the training and test sets which is completely forbidden!!!!!
            all_unique_training_labels = unique(cell2mat(the_training_label_numbers));
            all_unique_test_labels = unique(cell2mat(the_test_label_numbers));

    
            if length(intersect(all_unique_training_labels,  all_unique_test_labels)) > 0
                error('if running each split separately (i.e., use_unique_data_in_each_CV_split == 1), then none of the same labels can be in the training and test sets, otherwise there will be some of the same data in the training and test sets')
            end
        
            
            % remap labels (only using old 'test' data because I want each CV split to be independent)
            remapped_YTr_all = NaN .* ones(size(YTr_all));  % only getting labels (and data) from test set b/c want each CV split to contain unique data
            remapped_YTe_all = NaN .* ones(size(YTe_al));

            for iGroup = 1:length(the_training_label_numbers)
                remapped_YTr_all(ismember(YTr_all_cv, the_training_label_numbers{iGroup})) = iGroup;
                remapped_YTe_all(ismember(YTe_all, the_test_label_numbers{iGroup})) = iGroup;
            end
            
            
           % remove data from trials in which the labels are not appropriate for the training/test sets
           train_inds  = ~isnan(remapped_YTr_all);
           test_inds  = ~isnan(remapped_YTe_all);

           YTr_all = remapped_YTr_all(train_inds);   % remove NaNs from remapped labels
           YTe_all = remapped_YTe_all(test_inds);
                   
            for iCV = 1:length(XTr_all_time_cv{1})    
               for iTime = 1:length(XTr_all_time_cv) 
                   XTr_all_time_cv{iTime}{iCV} = XTr_all_time_cv{iTime}{iCV}(:, train_inds);   % only getting data from test set b/c want each CV split to be independent
                   XTe_all_time_cv{iTime}{iCV} = XTe_all_time_cv{iTime}{iCV}(:, test_inds);
               end
            end
   
           

        else
 
            
            % remap labels
            remapped_YTr_all = NaN .* ones(size(YTr_all));
            remapped_YTe_all = NaN .* ones(size(YTe_all));

            for iGroup = 1:length(the_training_label_numbers)
                remapped_YTr_all(ismember(YTr_all, the_training_label_numbers{iGroup})) = iGroup;
                remapped_YTe_all(ismember(YTe_all, the_test_label_numbers{iGroup})) = iGroup;
            end
            
            
           % remove data from trials in which the labels are not appropriate for the training/test sets
           training_data_inds_to_remove = isnan(remapped_YTr_all);
           test_data_inds_to_remove = isnan(remapped_YTe_all);

           for iTime = 1:length(XTr_all_time_cv)
               for iCV = 1:length(XTr_all_time_cv{1}) 
                    XTr_all_time_cv{iTime}{iCV}(:, training_data_inds_to_remove) = [];
                    XTe_all_time_cv{iTime}{iCV}(:, test_data_inds_to_remove) = [];  
               end
           end
           
           remapped_YTr_all(training_data_inds_to_remove) = [];
           remapped_YTe_all(test_data_inds_to_remove) = [];
           
           YTr_all = remapped_YTr_all;
           YTe_all = remapped_YTe_all;
             
        
        end
        
        
        
    end  % end get_data method
        
    
    
    
    
end   % end methods




end % end class
