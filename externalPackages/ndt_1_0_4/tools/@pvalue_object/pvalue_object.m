classdef pvalue_object < handle

%  This helper object calculates p-values from a file that has the decoding results in standard format, 
%   and a directory of files that has a number of 'null distribution' decoding results which were created by running 
%   the same decoding experiment a number of times with the labels randomly shuffled.  The p-values that are calculated
%   are based on a permutation test which gives the probability that the real decoding results came from 
%   the null distribution. The object also has a method to calculate the latency when the decoding results 
%   were first above chance (and it inherits from the handle class so it maintains its state after method calls).   
%
%  The constructor has the following form:
%   
%    pval_obj = pvalue_object(real_decoding_results_file_name, null_distribution_directory_name);
%
%    The input arguments to the constructor are: 
%
%       1. real_decoding_results_file_name:  A string specifying the name of a file that has real decoding results.
%           These results should be in standard results format (as created by the standard_resample_CV.run_cv_decoding method.
%
%       2.  null_distribution_directory_name: A string specifying the name of a directory that has multiple decoding 
%               result file that were run with the labels shuffled. The p-values are created by loading each file 
%               in this directory to create a null probability distribution that estimates what decoding accuracies would 
%               be expected to occur by chance.  
%
%     Both arguments to this constructor are optional, however the fields pval_obj.real_decoding_results_file_name and 
%       pval_obj.null_distribution_directory_name must both be set before calling the 
%       pval_obj.create_pvalues_from_nulldist_files method.    
%
%
%  The methods of this object are:  
%
%   1. [p_values null_distributions PVALUE_PARAMS] = pval_obj.create_pvalues_from_nulldist_files
%
%       This method creates p-values from the actual decoding results and the null distribution results (i.e., p-values are based on
%          a permutation test). Most specifically, all the files in the null distribution directory are loaded to create null distributions
%          at each point in time.  At each time point, a p-value is calculated based on the proportion of decoding values in the null distribution
%          that exceeds the real decoding result value. This method causes the properties pval_obj.p_values and pval_obj.null_distributions to
%          be set, and also returns the following output values:
%
%            p_values: a vector containing the p-values at each point in time (i.e., the probability one would get a decoding accuracy as high 
%               as the one reported if there was no relationship between the data and the class labels).
%
%            null_distributions: a [num_points_in_null_distribution x num_times] matrix containing the null distribution values at 
%               each point in time.
%
%            PVALUE_PARAMS: a list of parameters that were used to create the p-values (see the pval_obj.get_pvalue_parameters method below).
%
%
%   2.  [latency_time latency_ind] = pval_obj.get_latency_of_decoded_information 
%
%       This method returns an estimate of the latency (i.e., time) when the decoding results are above chance. The  
%           pval_obj.create_pvalues_from_nulldist_files method must be called (or the field  pval_obj.p_values must
%           be manually set) prior to running this method, and additionally latency_time_interval must also be set.
%           This method works by finding all p_values that are less than or equal to pval_obj.latency_alpha_significance_level, 
%           and then returning the first time in which the p-values are at or below this significance level for 
%           pval_obj.latency_num_consecutive_sig_bins.
%           
%
%   3.  PVALUE_PARAMS = pval_obj.get_pvalue_parameters
%
%       This method returns a number of parameters that were used to calculate the p-values and get the decoding information latency. The
%           structure, PVALUE_PARAMS, that is returned by this object has the following fields:
%
%               .num_points_in_null_distribution:  The number of points in the null distribution.  
%               
%               .smallest_significance_level:  the smallest possible significance level that can be used based on the number of points 
%                   in the null distribution (i.e., this value is equal to 1/PVALUE_PARAMS.num_points_in_null_distribution).  If the 
%                   parameter pval_obj.latency_alpha_significance_level is set to 0 (default value), then the p-values should be listed as 
%                   p < PVALUE_PARAMS.smallest_significance_level.  
%               
%               .result_type_name: A string specifying what type of results were used to create these p-values (e.g., ZERO_ONE_LOSS, etc.)
%
%               .training_time_ind_to_use:  The pval_obj.training_time_ind_to_use value that was used.
%
%               .real_decoding_results_file_name:  The pval_obj.real_decoding_results_file_name string that was used.  
%        
%               .null_distribution_directory_name:  The pval_obj.null_distribution_directory_name string that was used.
%
%               .saved_results_structure_name:  The pval_obj.saved_results_structure_name string that was used.
%
%               .null_distribution_file_prefix_name: The pval_obj.null_distribution_file_prefix_name that was used.
%
%               .latency_alpha_significance_level: The pval_obj.latency_alpha_significance_level value that was used.
%
%               .latency_num_consecutive_sig_bins = The pval_obj.latency_num_consecutive_sig_bins that was used.
%
%               .latency_time_interval = The pval_obj.latency_time_interval that was used.
%
%
%
%   The following properties can be set to change the behavior of this object:
%
%       1. real_decoding_results_file_name: the name of the non-shuffled decoding results that are compared to the
%           null distribution to see when the results are above chance. This value can also be set in the constructor.
%
%       2. null_distribution_directory_name: the name of the directory that contains the null distribution results files.
%           This value can also be set in the constructor.    
%
%       3. the_result_type: (default = 1).  Specifies which result type should be used to create p-values from.
%           If this is set to 1, then 0-1 loss results are plotted.
%           If this is set to 2, normalized rank results are plotted. 
%           If this is set to 3, the mean decision values are plotted. 
%           If this is set to 4, ROC AUC results run separately on each CV split are plotted.    
%           If this is set to 5, ROC AUC results combined over CV splits are plotted.    
%           If this is set to 6, mutual information created from a confusion matrix that combining data from all resamples is plotted.
%           If this is set to 7, mutual information created from a confusion matrix that is calculated separate and then averaged over resamples is plotted.    
%
%       4. null_distribution_file_prefix_name (deafult is '').  A string that specifies what the beginning of the names of files in the 
%           null distribution directory is.  This is useful if there are multiple types of results (e.g., from different decoding analyses)
%           stored in the directory that has the null files but you only want the results in some of these files to be used.  
%
%       5. training_time_ind_to_use (default = -1).  If a full TCT plot has been created, this specifies which row (i.e., training time bin) 
%           of the TCT plot should be used when calculating the p-values.  Setting this to a value of less than zero creates p-values when the classifier 
%           was trained and tested from the same time bin (i.e., the diagonal of the TCT plot, or equivalent vector of results if the classifier was only trained and 
%           tested at the same time).  
%
%       6. saved_results_structure_name (default is 'DECODING_RESULTS').  A string specifying the name of the variable that has the decoding results.
%
%       7. p_values (default is []).  These are the p-values that are usually set by calling the pval_obj.create_pvalues_from_nulldist_files method.
%           However, one can set these values manually, and then use the pval_obj.get_latency_of_decoded_information method.  Doing this is useful if
%           one is getting the latency many times so that one calculate the p-values once, save them, and then load them into this object to get the latency.
%
%       8. latency_time_interval:  This specifies which time points in the experiment the p-values correspond to, and must be set prior to calling the 
%           pval_obj.get_latency_of_decoded_information method.  This property can be set to either: 1) a vector specifying which times to use, 
%           2) a time_interval_object that one can get time_intervals from or 3) structure with the fields latency_time_interval.bin_width, 
%           latency_time_interval.sampling_interval, and optionally latency_time_interval.alignment_event_time which will create a time interval
%           that has the corresponding bin widths, step sizes and zero time, and latency_time_interval.start_time and latency_time_interval.end_time
%           which will give the start and end times of the time interval. 
%
%       9. latency_time_interval_alignment_type: (default = [], which means show the beginning and end of the first significant time bin). If latency_time_interval
%            is set to a vector of numbers, this property will be ignored.  However if latency_time_interval is set to a time_interval_object or to a structure
%            containing latency_time_interval.bin_width and latency_time_interval.sampling_interval then this will cause the latency estimate to report either: 
%            the beginning time of the first significant time bin (latency_time_interval_alignment_type = 1), the middle time of the first significant time bin 
%            (latency_time_interval_alignment_type = 2), the end time of the first significant time bin (latency_time_interval_alignment_type = 3), the beginning 
%            and end time of the first significant time bin (latency_time_interval_alignment_type = 4, default), the time interval of data that was added to 
%            make a bin significant, relative to the previous bin which was not significant (latency_time_interval_alignment_type = 5), or use the alignment 
%            already specified in the time_interval_object or by the structure latency_time_interval.alignment time (latency_time_interval_alignment_type = 6).
%
%      10. latency_alpha_significance_level (default = 0).  The significance level (alpha value) that the p-values must be less than to be considered to have not 
%           occurred by chance. When calculated the latency of decoding information, all the p-values are compared to this significance level value to determine
%           whether a time point is considered significant (and the latency is determined based on whether there are pval_obj.latency_num_consecutive_sig_bins
%           significant bins in a row).  
%
%      11. latency_num_consecutive_sig_bins (default = 5).  The number of consecutive time bins that must be significant in order for a specific time to be
%           considered the time point when the results are first above chance.  The reason this property is needed is because calculating p-values at many
%           time periods introduces a high probability that one will have a small p-value even when the null hypothesis is correct (this is due to multiple
%           comparisons issues that commonly effect null-hypothesis signifiance tests).  A commonly used (ad hoc) method in neuroscience to deal 
%           with this issue is to define the latency as the time when the p-values remain significant of multiple consecutive bins, which is the method we 
%           are using here (empirically it appears to produce reasonable results).  
%
%       12. real_decoding_results_lower_than_null_distribution (default = 0).  If this is set to one, the p-values are calculated based on the proportion of 
%            null distribution decoding results are lower than the actual real decoding result (i.e., the test shows the probability that the real decoding
%            result would have been that low by chance).  
%
%       13. collapse_all_times_when_estimating_pvals (default = 0). If this is set to one, the null distributions from all time bins are combined together to create 
%            one larger total null distribution.  The p-values are then calculated by comparing the actual decoding accuracy at each point 
%            in time to this larger null distribution (with this same null distribution is used for all points in time). The advantage of using this is that if 
%            the null distributions at each point in time are the same, then one can get a more precise estimate of the p-values for the same computational cost. [Added with NDT version 1.4] 
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
      
      real_decoding_results_file_name = [];  % name of the file that has the real decoding results
      null_distribution_directory_name = []; % name of the directory that has the null distribution results files   
      the_result_type = 1;                  % the type of decoding results to use to get the p-values (default, use 0-1 loss results)
      null_distribution_file_prefix_name = '';  % prefix string specifying to only use null distribution files that start with this string
      training_time_ind_to_use = -1;  % if this is less than 1, the diagonal of TCT matrix is used (or if mean decoding results is a vector, then this vector is used); otherwise this ind in the TCT plot is used.
      saved_results_structure_name = 'DECODING_RESULTS';  % name of variable that has the saved standard results
    
      p_values = [];  % the p-values.  This can be set using the create_pvalues_from_nulldist_files, or can set it directly to a vector and then use the get_latency method to calculate information latency
      
      
      % information needed to an estimate of the latency when the results are above chance
      latency_alpha_significance_level = 0;  % alpha level that the p-values must be less than to be considered significant (i.e., decoding results above chance).
      latency_num_consecutive_sig_bins = 5;  % the number of consecutive time bins that the results must be significant in order for a time period to be considered the correct information latency.
      latency_time_interval;  % the time intervals used for calculating the latency of information        
      
      latency_time_interval_alignment_type = [];  % if the latency_time_interval is set to a binning parameters structure, or to a time_interval object, then this parameter can be set to return different estimates of the latency (see documentaiton above)
                                   % if empty, a value of 4 is the default 
     
      real_decoding_results_lower_than_null_distribution = 0;  % if set to 1, will calculate p-value based on proportion of null distribution results less than real decoding result                             
                                   
      collapse_all_times_when_estimating_pvals = 0;  % if set to 1, will collapse all the null distribution values from all points in time into a single null distribution.
                                                        % this larger null distribution will then be used when calculating p-values at all points in time. 
  end
    
  
  
  properties (GetAccess = 'public', SetAccess = 'private')
      null_distributions;  % the null distributions calculated from loading all the null distribution files
      real_results;        % the real decoding results calculated from loading the real decoding result file
      num_points_in_null_distribution;  % the number of points in the null distribution (which is the same as the number of files in the null-distribution directory)
      result_type_name; %  the name of the type of results are that are used to calculate the p-values (i.e., ZERO_ONE_LOSS_RESULTS, etc.).
  end
  
  
  
  
  
  
  methods

  
        % the constructor (two arguments can be set here, or they can be set later prior to calling the create_pvalues_from_nulldist_files method)     
        function pval_obj = pvalue_object(real_decoding_results_file_name, null_distribution_directory_name)
            
            if nargin > 0
                pval_obj.real_decoding_results_file_name = real_decoding_results_file_name;
            end
            
            if nargin > 1
                pval_obj.null_distribution_directory_name = null_distribution_directory_name; 
            end
            
        end
        
        
        
        
        
        % the main method that creates the p-values from a bunch of null distribution files
        function [p_values null_distributions PVALUE_PARAMS] = create_pvalues_from_nulldist_files(pval_obj)

            
            % sanity checks to make sure proper fields are set before this method is called
            if isempty(pval_obj.real_decoding_results_file_name)
                error('pval_obj.real_decoding_results_file_name must be set to a string specifying a file that has standard results prior to calling the pval_obj.create_pvalues_from_nulldist_files method');
            end
                       
            if isempty(pval_obj.null_distribution_directory_name)
                error('pval_obj.null_distribution_directory_name must be set to a string specifying a directory that has null distribution decoding files prior to calling the pval_obj.create_pvalues_from_nulldist_files method');
            end
            
            
            
            % determine which type of decoding results should be used
            if pval_obj.the_result_type == 1
                 pval_obj.result_type_name = 'ZERO_ONE_LOSS_RESULTS';
            elseif pval_obj.the_result_type == 2
                 pval_obj.result_type_name = 'NORMALIZED_RANK_RESULTS';
            elseif pval_obj.the_result_type == 3
                 pval_obj.result_type_name = 'DECISION_VALUES';
            elseif pval_obj.the_result_type == 4
                 pval_obj.result_type_name = 'ROC_AUC_RESULTS.separate_CV_ROC_results';
            elseif pval_obj.the_result_type == 5
                 pval_obj.result_type_name = 'ROC_AUC_RESULTS.combined_CV_ROC_results';
            elseif pval_obj.the_result_type == 6    
                 pval_obj.result_type_name = 'MUTUAL_INFORMATION.from_combined_confusion_matrix_over_all_resamples';
            elseif pval_obj.the_result_type == 7    
                 pval_obj.result_type_name = 'MUTUAL_INFORMATION.from_separate_confusion_matrix_for_each_resample';                  
            end

            
            
            % load the real decoding results
            all_real_data = load(pval_obj.real_decoding_results_file_name);
            
            if pval_obj.the_result_type == 6
                pval_obj.real_results = eval(['all_real_data.' pval_obj.saved_results_structure_name '.' pval_obj.result_type_name '.decoding_results']);
            else
                pval_obj.real_results = eval(['all_real_data.' pval_obj.saved_results_structure_name '.' pval_obj.result_type_name '.mean_decoding_results']);
            end
            

            % specify which row of the TCT matrix to use so that p-values for a particular training time period can potentially be returned
            if pval_obj.training_time_ind_to_use < 1 
                if size(pval_obj.real_results, 2) > 1
                    pval_obj.real_results = diag(pval_obj.real_results); 
                end
            else
               pval_obj.real_results = real_data(pval_obj.training_time_ind_to_use, :);
            end



            % create the null distribution for each time bin
            the_null_dist_dir = dir([pval_obj.null_distribution_directory_name pval_obj.null_distribution_file_prefix_name '*.mat']);

            for iNullDist = 1:length(the_null_dist_dir)

                curr_all_null_data = load([pval_obj.null_distribution_directory_name the_null_dist_dir(iNullDist).name]);   
                
                if pval_obj.the_result_type == 6
                   curr_null_results = eval(['curr_all_null_data.' pval_obj.saved_results_structure_name '.' pval_obj.result_type_name '.decoding_results']);
                else
                   curr_null_results = eval(['curr_all_null_data.' pval_obj.saved_results_structure_name '.' pval_obj.result_type_name '.mean_decoding_results']);
                end

                if pval_obj.training_time_ind_to_use < 1 
                    if size(curr_null_results, 2) > 1
                        curr_null_results = diag(curr_null_results ); 
                    end
                else
                   curr_null_results = curr_null_results(pval_obj.training_time_ind_to_use, :);
                end

                pval_obj.null_distributions(iNullDist, :) = curr_null_results;

            end

            
            
            if pval_obj.collapse_all_times_when_estimating_pvals == 1       % added with NDT version 1.4 - creates null distribution using data points from all points in time
                
                total_null_distribution_collapsed_across_time = pval_obj.null_distributions(:);
                              
                for iTime = 1:size(pval_obj.null_distributions, 2)

                    if pval_obj.real_decoding_results_lower_than_null_distribution   % if real_decoding_results_lower_than_null_distribution == 1, see if real results are less than null distirbution values
                        pval_obj.p_values(iTime) = length(find(total_null_distribution_collapsed_across_time <= pval_obj.real_results(iTime)))./length(total_null_distribution_collapsed_across_time);   
                    else
                        pval_obj.p_values(iTime) = length(find(total_null_distribution_collapsed_across_time >= pval_obj.real_results(iTime)))./length(total_null_distribution_collapsed_across_time);
                    end
                end
            
                pval_obj.num_points_in_null_distribution = numel(total_null_distribution_collapsed_across_time);  % how many points are in the null distribution 
                
            else
                
                % create p-values from the null_distribution at each time point
                for iTime = 1:size(pval_obj.null_distributions, 2)

                    if pval_obj.real_decoding_results_lower_than_null_distribution   % if real_decoding_results_lower_than_null_distribution == 1, see if real results are less than null distirbution values
                        pval_obj.p_values(iTime) = length(find(pval_obj.null_distributions(:, iTime) <= pval_obj.real_results(iTime)))./size(pval_obj.null_distributions, 1);   
                    else
                        pval_obj.p_values(iTime) = length(find(pval_obj.null_distributions(:, iTime) >= pval_obj.real_results(iTime)))./size(pval_obj.null_distributions, 1);   
                    end

                end

                pval_obj.num_points_in_null_distribution = size(pval_obj.null_distributions, 1);  % how many points are in the null distribution 
                
            end
                
  
            % return the results
            p_values = pval_obj.p_values;
            null_distributions = pval_obj.null_distributions;
            PVALUE_PARAMS = get_pvalue_parameters(pval_obj);
            
            
        end

        
        
        % a function that returns parameters used to create the null-distribution and get the information latency
        function PVALUE_PARAMS = get_pvalue_parameters(pval_obj)
            
            
            % parameters used to generate the p-values from the null distribution
            PVALUE_PARAMS.num_points_in_null_distribution = pval_obj.num_points_in_null_distribution;
            PVALUE_PARAMS.smallest_significance_level = 1./pval_obj.num_points_in_null_distribution;
            PVALUE_PARAMS.result_type_name = pval_obj.result_type_name;
            PVALUE_PARAMS.training_time_ind_to_use = pval_obj.training_time_ind_to_use;
            
            PVALUE_PARAMS.real_decoding_results_file_name = pval_obj.real_decoding_results_file_name;
            PVALUE_PARAMS.null_distribution_directory_name = pval_obj.null_distribution_directory_name; 
            PVALUE_PARAMS.saved_results_structure_name = pval_obj.saved_results_structure_name;
            PVALUE_PARAMS.null_distribution_file_prefix_name = pval_obj.null_distribution_file_prefix_name;
            
            
            % parameters for getting the latency of when information was above chance
            PVALUE_PARAMS.latency_alpha_significance_level = pval_obj.latency_alpha_significance_level;
            PVALUE_PARAMS.latency_num_consecutive_sig_bins = pval_obj.latency_num_consecutive_sig_bins;
            PVALUE_PARAMS.latency_time_interval = pval_obj.latency_time_interval;           

        end
        
        
        
        
        % a function that returns the latency when decoded information was first above chance
        function [latency_time latency_ind] = get_latency_of_decoded_information(pval_obj) 
            
            
            % find the bin index where the p-values are first below chance
            
            % definition of latency:  the first time p-values < the alpha significance level, and stays below that level for latency_num_consecutive_sig_bins time periods
            nonsignificant_pvalue_times = double(~(pval_obj.p_values <= pval_obj.latency_alpha_significance_level));  % 1 indicates nonsignicant time, 0 indicates significant time
           
            conv_consec_sig_times = conv(ones(pval_obj.latency_num_consecutive_sig_bins, 1), nonsignificant_pvalue_times);
            conv_consec_sig_times = conv_consec_sig_times(pval_obj.latency_num_consecutive_sig_bins:end);   
                        
            latency_ind = min(find(conv_consec_sig_times == 0));
            
            % if there is no time period that is statistically significant, return NaN
            if isempty(latency_ind)
                latency_time = NaN;
                return
            end
            
            
            % calculate what time that corresponds to (which depends on the time interval definition given)

            
            % if time_interval_object or a structure that one can create a time interval object, create the time interval object from it            
           if (isobject(pval_obj.latency_time_interval) && strcmp(class(pval_obj.latency_time_interval), 'time_interval_object')) || isstruct(pval_obj.latency_time_interval)  

               
              if isstruct(pval_obj.latency_time_interval) 
                time_interval_obj = time_interval_object;
                time_interval_obj.set_parameters_with_binning_stucture(pval_obj.latency_time_interval);
              else
                time_interval_obj = pval_obj.latency_time_interval;
              end
              
              latency_time_interval_alignment_type = pval_obj.latency_time_interval_alignment_type;
              if isempty(latency_time_interval_alignment_type), latency_time_interval_alignment_type = 4; end
              
              
              if (latency_time_interval_alignment_type > 0) && (latency_time_interval_alignment_type < 4)   % algin to begining (1), middle (2), or end (3) of the intervals
                  
                  time_interval_obj.alignment_type = latency_time_interval_alignment_type;
                  curr_latency_time_intervals = time_interval_obj.get_time_interval;
                  latency_time = curr_latency_time_intervals(latency_ind);
              
              elseif latency_time_interval_alignment_type == 4  % return both the beginning and the end of the interval
              
                  time_interval_obj.alignment_type = 1;
                  curr_latency_time_intervals = time_interval_obj.get_time_interval;
                  latency_time(1) = curr_latency_time_intervals(latency_ind);
                  
                  time_interval_obj.alignment_type = 3;
                  curr_latency_time_intervals = time_interval_obj.get_time_interval;
                  latency_time(2) = curr_latency_time_intervals(latency_ind);

              elseif latency_time_interval_alignment_type == 5   % return the amount of time added that made the results fall below chance
                   
                  time_interval_obj.alignment_type = 3;
                  curr_latency_time_intervals = time_interval_obj.get_time_interval;
                  latency_time(2) = curr_latency_time_intervals(latency_ind);
                  latency_time(1) = latency_time(2) -  time_interval_obj.sampling_interval;  
                  
              elseif latency_time_interval_alignment_type == 6   % use the alignment value given in the binning_parameters structure or in the time_interval_object

                  curr_latency_time_intervals = time_interval_obj.get_time_interval;
                  latency_time = curr_latency_time_intervals(latency_ind);
                  
              else
                    error('The property latency_time_interval_alignment_type must be set to a value from 1 to 6 or be empty')
              end
                  
              
              curr_latency_time_interval = time_interval_obj.get_time_interval;
                             

           elseif isvector(pval_obj.latency_time_interval)
                curr_latency_time_intervals = pval_obj.latency_time_interval;
                latency_time = curr_latency_time_intervals(latency_ind);
           else
                error('pval_obj.latency_time_interval must be set to either a vector of times, or a structure that has the fields bin_width and sampling_interval before calling get_latency_of_decoded_information')
           end
            
                       
        end
        

        
        
        
    end



end



























