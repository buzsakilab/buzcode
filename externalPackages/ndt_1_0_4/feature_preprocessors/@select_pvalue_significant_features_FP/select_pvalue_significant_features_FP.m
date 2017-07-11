classdef select_pvalue_significant_features_FP
    
% select_pvalue_significant_features_FP is a feature preprocessing (FP) object 
% that uses an ANOVA to find the features with the lowest p-values using data
% from the training set.  The algorithm uses all features that are less than a 
% set p-value threshold.
%
%
% Like all FP objects, there are three main methods, which are:
%
%  1.  [fp XTr_preprocessed] = set_properties_with_training_data(fp, XTr, YTr) 
%        This method takes the training data (XTr) and calculates the most selective
%        features using an ANOVA (i.e., features with the smallest pvalues), and it
%        changes the state of the object to list which features were selected.  It 
%        then returns a version of the training data (XTr_fp_processed), that only
%        has features that have p-values less than a set threshold.
% 
%  2.  X_preprocessed = preprocess_test_data(fp, X_data)
%         This method takes other data (e.g., the test data), and returns the data
%         that has only the features pvalue found using the training set.
%
%  3.  current_FP_info_to_save = get_current_info_to_save(fp).  Returns
%       the p-values from ANOVAs applied to each feature/neuron in the structure
%       current_preprocessing_information_to_save.the_p_values_org_order that will
%       be saved the CV algorithm.
%
%
% Additionally, this object has the following properties that can be set
%  
%  1. pvalue_threshold
%       Select the alpha level that determined whether a site is considered significant and should be 
%       included as one of the features that is used.  This must be set prior to running this code.
%      
%  2. save_extra_info 
%       Causes the get_current_info_to_save method to return the p-values of each feature
%       which will be saved the CV algorithm. 
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
        pvalue_threshold = [];
        save_extra_info = 0;  % flag if one wants to save extra information about feature pvalues
    end

          
    properties  (GetAccess = 'public', SetAccess = 'private')
        feature_inds_to_use = [];
    end
    
    
    properties (Access = private)
        current_preprocessing_information_to_save = [];
    end
       
    
    
    
    methods 

        % constructor 
        function fp = select_pvalue_significant_features_FP
        end

        
        function current_FP_info_to_save = get_current_info_to_save(fp)
            current_FP_info_to_save = fp.current_preprocessing_information_to_save;
        end
        

        
        function [fp XTr_preprocessed] = set_properties_with_training_data(fp, XTr, YTr)   
            
            
            if isempty(fp.pvalue_threshold)  % sanity check 
                error('pvalue_threshold must be set prior to calling this method')
            end
            
            
            [sorted_features ANOVA_pvals] = rank_features_using_an_ANOVA(XTr, YTr);   % much faster way to do a one-way ANOVA on all features
            
            % these parameters are not used by the algorithm, but they might be interesting to save anyway... 
            if fp.save_extra_info == 1
                %fp.current_preprocessing_information_to_save.neuron_rank = sorted_features;  % don't really need to save these, can sort pvalues later if I want this info
                %fp.current_preprocessing_information_to_save.the_p_values_sorted = sorted_p_values;  
                fp.current_preprocessing_information_to_save.the_p_values_org_order = ANOVA_pvals;
            end
            

            fp.feature_inds_to_use = find(ANOVA_pvals < fp.pvalue_threshold);
            
            % when none of the neurons are less than the specified threshold, only use the single most selective neurons
            %  (this is useful so that the algorithm does not crash during the baseline period before any of the neurons should be selective)
            if length(fp.feature_inds_to_use) < 1    
                warning('there are no significant features - using just one neuron that has the smallest p-value!!!!!!!')
                [val fp.feature_inds_to_use] = min(ANOVA_pvals);   % should always be the first ind???
            end
            
                        
            XTr_preprocessed = XTr(fp.feature_inds_to_use, :);
               
        end

        

        function X_preprocessed = preprocess_test_data(fp, X_data)
            X_preprocessed = X_data(fp.feature_inds_to_use, :); 
        end


        
    end  % end methods
    
end    % end classdef




    
    
    













