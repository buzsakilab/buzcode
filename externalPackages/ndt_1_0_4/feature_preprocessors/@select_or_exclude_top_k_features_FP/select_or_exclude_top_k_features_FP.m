classdef select_or_exclude_top_k_features_FP
    
% select_or_exclude_top_k_features_FP is a feature preprocessing (FP) object 
%   that uses an ANOVA to find the k features with the lowest p-values using data
%   from the training set.  The algorithm either uses only these features 
%   (i.e., removes all other features from the training and test sets), or 
%   removes these features (i.e., returns training and test sets that have 
%   these features removed).  The object can also exclude the top k features, and then
%   use the next best j features.  
%   
%
%
% Like all FP objects, there are three main methods, which are:
%
%  1. [fp XTr_preprocessed] = set_properties_with_training_data(fp, XTr, YTr) 
%       This method takes the training data (XTr) and calculates the most selective
%       features using an ANOVA (i.e., features with the smallest pvalues), and it
%       changes the state of the object to list which features were selected.  It 
%       then returns a version of the training data (XTr_fp_processed), that either
%       only have these features, or have these features removed.
% 
%  2. X_preprocessed = preprocess_test_data(fp, X_data)
%       This method takes other data (e.g., the test data), and either uses
%       or excludes these features from the data.  
%
%  3. current_FP_info_to_save = get_current_info_to_save(fp).  Returns
%       the p-values from ANOVAs applied to each feature/neuron in the structure
%       current_preprocessing_information_to_save.the_p_values_org_order that will
%       be saved the CV algorithm.
%
%
% Additionally, this object has the properties that can be set
%  
%  1. num_features_to_exclude
%       Select how many of the most selective features to exclude.
%
%  2. num_features_to_use
%       Select how many of the most selective features to use (will exclude all other features less selective).
%       Either num_features_to_exclude and/or num_features_to_use must be set for this feature-preprocessor to work.
%      
%  3. save_extra_info
%      Causes the get_current_info_to_save method to return the p-values of each feature
%      which will be saved the CV algorithm.  
%
%  Note: Both num_features_to_exclude and num_features_to_use can be set to non-zero values.  
%  If this is done then the most selective num_features_to_exclude will be excluded
%  and then the next num_features_to_use selective features will be used.  
%  This allows one to do things like compare the classification accuracy using
%  only the top k features, to the classification accuracy for using the next
%  k most selective features, etc. (e.g., to look at the drop-off in classification
%  accuracy as using progressively less selective features.
%
%  Note: X_data is in the form [num_features x num_examples] 
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
        num_features_to_exclude = [];
        num_features_to_use = [];
        save_extra_info = 0;  % flag if one wants to save extra information about feature pvalues
    end

    
    properties (GetAccess = 'public', SetAccess = 'private')
        feature_inds_to_use = [];  % the main variable that is set by this FP algorithm which indicates which features should be used
    end
    
    
    properties (Access = private)
        current_preprocessing_information_to_save = [];
    end
    
    
    
    methods 
    
        % constructor 
        function fp = select_or_exclude_top_k_features_FP
        end

        
        function current_FP_info_to_save = get_current_info_to_save(fp)
            current_FP_info_to_save = fp.current_preprocessing_information_to_save;
        end
        

        function [fp XTr_preprocessed] = set_properties_with_training_data(fp, XTr, YTr)   
            
            
            [sorted_features ANOVA_pvals] = rank_features_using_an_ANOVA(XTr, YTr);   % much fast way to do a one-way ANOVA on all features

            
            % these parameters are not used by the algorithm, but they might be interesting to save anyway... 
            if fp.save_extra_info == 1
                %fp.current_preprocessing_information_to_save.neuron_rank = sorted_features;  % don't really need to save these, can sort pvalues later if I want this info
                %fp.current_preprocessing_information_to_save.the_p_values_sorted = sorted_p_values;  
                fp.current_preprocessing_information_to_save.the_p_values_org_order = ANOVA_pvals;
            elseif fp.save_extra_info ~= 0
                error('fp.save_extra_info must be set to either 0 or 1')
            end
            

            fp.feature_inds_to_use = sorted_features;
            
            if (isempty(fp.num_features_to_exclude)) && (isempty(fp.num_features_to_use))    % sanity check to make sure some features are used/excluded, otherwise there is no point in using this object
                error('either num_features_to_exclude or num_features_to_use must be set prior to calling this function')
            end
            
            
            % if both fp.num_features_to_exclude and fp.num_features_to_use are set, then first exclude the best features and then use the next best features...
                       
            if ~isempty(fp.num_features_to_exclude)
                fp.feature_inds_to_use = fp.feature_inds_to_use((fp.num_features_to_exclude + 1):end);
            end
            
            if ~isempty(fp.num_features_to_use)
               fp.feature_inds_to_use = fp.feature_inds_to_use(1:fp.num_features_to_use);
            end
            
            
            XTr_preprocessed = XTr(fp.feature_inds_to_use, :);
               
                        
        end

        

        function X_preprocessed = preprocess_test_data(fp, X_data)
            X_preprocessed = X_data(fp.feature_inds_to_use, :); 
        end


        
    end  % end methods
    
end    % end classdef




    
    
    













