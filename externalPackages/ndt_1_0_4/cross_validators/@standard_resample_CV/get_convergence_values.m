function [maximum_absolute_difference_between_mean_and_ith_deleted_result percentage_relative_to_max_between_mean_and_ith_result_deleted] = get_convergence_values(cv, results_for_each_resample_run_averaged_over_CV_splits)


 % This code calculates the different between the mean results averaged over all resample runs, and the mean result with the ith resample run deleted 
 %   (for all i resample runs). It then find the maximum difference between this overall mean and ith resample deleted mean, with the maximum taken over all
 %   resample runs, and training and test times (full TCT matrix). If this maximum difference is less than the threshold specified in 
 %   cv.stop_resample_runs_when_decoding_results_change_by_less_than then we return that the results have converged (results_have_converged = 1) and results_have_converged = 0
 %   otherwise.  The maximum difference value between the mean results and the ith resample deleted mean results in the variable 
 %   maximum_different_between_mean_and_ith_deleted_result. Put more simply, this function calculates the largest effect that any single resample run had on the mean
 %   decoding results at any point in time. 
 %
 %  Outputs: 
 %    1. maximum_absolute_difference_between_mean_and_ith_deleted_result: the maximum value that the mean of a decoding measure changed when the ith resample run is deleted
 %       (where the maximum is taken over all resample runs and time points)
 %
 %    2. percentage_relative_to_max_between_mean_and_ith_result_deleted: shows how much deleting the most influential single resample run affects the mean results, relative to the best
 %        mean decoding performance. Since plots are typically scaled to the maximum decoding value, this shows the percentage of variance there is in the plot that could potentially 
 %        be reduced if more resampling runs were used (i.e., would the plot look less wiggly if more resample runs were used). 
 
 
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


num_resample_runs_completed = size(results_for_each_resample_run_averaged_over_CV_splits, 1);

% sum of results over all resample runs (is equivalent to the averaged results * num_resample_runs_completed)
sum_results_over_resample_runs = squeeze(sum(results_for_each_resample_run_averaged_over_CV_splits, 1));   

% created matrix where summed results are repeated num_resample_runs_completed times, i.e., result is of size [num_resample_runs_completed num_train_times num_test_times]
repeated_summed_result = squeeze(shiftdim(repmat(sum_results_over_resample_runs, [1, 1, num_resample_runs_completed]), 2));

% find the mean of all the results when the results from the ith resample run is deleted
mean_results_ith_resample_run_deleted = (repeated_summed_result - results_for_each_resample_run_averaged_over_CV_splits)./(num_resample_runs_completed - 1);


% see how much the overall mean differents from the mean calculated when the ith resample run is deleted
overall_mean_results_minus_ith_resample_deleted_result = (repeated_summed_result./num_resample_runs_completed) - mean_results_ith_resample_run_deleted;


% find the maximum absolute difference between the mean results and the mean results when the ith resample result is deleted -> over all times and resample runs 
maximum_absolute_difference_between_mean_and_ith_deleted_result = max(max(max(abs(overall_mean_results_minus_ith_resample_deleted_result))));


% find the maximum percentage difference between the mean results and the mean results when the ith resample result is deleted -> over all times and resample runs 
% this is a poor measure because some types of decoding results converage to zero during baseline (e.g., decision values with the MCC classifier) making this statistic skewed...
% maximum_percentage_change_between_mean_and_ith_deleted_result = max(max(max(abs(overall_mean_results_minus_ith_resample_deleted_result./(repeated_summed_result./num_resample_runs_completed))))) .* 100; % (over_mean - ith_deleted)./overall_mean;


% shows how much the largest deviation with the ith resampole is deleted differs from the max of the (mean) decoding accuracy over 
% (since plot will typically be scaled to max, will show how large the deviations are relative to the scale of the plot)
percentage_relative_to_max_between_mean_and_ith_result_deleted = (maximum_absolute_difference_between_mean_and_ith_deleted_result./max(max((mean(results_for_each_resample_run_averaged_over_CV_splits, 1))))) .* 100;













% 
% 
% 
% 
% 
% 
% 
% % junk
% 
% 
% 
% 
% i = 1;
% num_training_times = 0;
% while num_training_times(1) == 0
%     if i ~= 4
%         num_training_times = size(results_to_display{i}, 3);
%         num_test_times = size(results_to_display{i},4);
%     else
%         num_training_times = size(results_to_display{i}, 4);
%         num_test_times = size(results_to_display{i}, 5);
%     end
%     i = i + 1;
% end
% 
% 
% 
%  % default behavior if there is only 1 training point and display_progress.training_time_to_display_results == -1, is to just use this one training time for the results
%  if (display_progress.training_time_to_display_results == -1)  && (num_training_times == 1)
%      display_progress.training_time_to_display_results = 1;
%  end
% 
% 
%  % display_progress.training_time_to_display_results should be set if the num_training_times ~= num_test_times
%  if (display_progress.training_time_to_display_results == -1) && (num_training_times~= num_test_times) && (num_test_times ~= 1)
% 
%      warning(['If the number of training times does not equal the number of test times, and one wants to display the resampled CV progress, ' ...
%                 'then display_progress.training_time_to_display_results should be set.  Using the default value of display_progress.training_time_to_display_results = 1, ' ...
%                 'to display the results here.'])
% 
%      display_progress.training_time_to_display_results = 1;
% 
%  end
%  
%  
%  
%  
% current_results = [];
%  
% 
%  for iResult = 1:length(results_to_display)
%        
%               
%      iResample = size(results_to_display{1}, 1);
% 
%      if iResample == 0
%          continue;
%      end
%      
% 
%     % if a particular training time that has been specified to display the results
%     % note: if only saving results when training and testing at the same time period, then the displayed results will only be for one time point 
%     if display_progress.training_time_to_display_results > 0
%           
% 
%         if iResult ~= 4
%            curr_bootstrap_progress_results = squeeze(mean(results_to_display{iResult}(:, :, display_progress.training_time_to_display_results, :), 2));  
%         else
%            curr_bootstrap_progress_results = squeeze(mean(mean(results_to_display{iResult}(:, :, :, display_progress.training_time_to_display_results, :), 3), 2));  
%         end
%         
%         if iResample == 1
%            curr_bootstrap_progress_results = curr_bootstrap_progress_results';
%         end
%         
%     else
%         
%         % if only saving the results for training and testing at the same time
%         if ((length(size(results_to_display{iResult})) == 3) && (iResult ~=4)) || ((length(size(results_to_display{iResult})) == 4)  && (iResult == 4))
%             
%             if iResult ~=4
%                  curr_bootstrap_progress_results = squeeze(mean(results_to_display{iResult}, 2));
%             else
%                  curr_bootstrap_progress_results = squeeze(mean(mean(results_to_display{iResult}, 3), 2));
%             end
%             
%             if iResample == 1
%                curr_bootstrap_progress_results = curr_bootstrap_progress_results';
%             end            
%             
%         else
%         
%              for iDisplayResultTrainingTime = 1:num_training_times 
% 
%                  if iResult ~= 4
%                     curr_bootstrap_progress_results(:, iDisplayResultTrainingTime) = squeeze(mean(results_to_display{iResult}(:, :, iDisplayResultTrainingTime, iDisplayResultTrainingTime), 2));
%                  else
%                     curr_bootstrap_progress_results(:, iDisplayResultTrainingTime)  = squeeze(mean(mean(results_to_display{iResult}(:, :, :, iDisplayResultTrainingTime, iDisplayResultTrainingTime), 3), 2));
%                  end
%              end
%          
%         end
%          
%          
%          
%     end
%         
%         
%      
%    if iResample == 1
%        current_results = [current_results; curr_bootstrap_progress_results; NaN .* ones(size(curr_bootstrap_progress_results))];   % can't show cofficients of variation if only 1 bootstrap has been run (so skip the following part)
%    else
% 
%        curr_total_mean = mean(curr_bootstrap_progress_results, 1);
%        curr_one_bootstrap_left_out_mean = ((repmat(curr_total_mean, [iResample 1]) - (curr_bootstrap_progress_results./iResample)) .* (iResample./(iResample -1)))';   % means over bootstraps with the ith bootstrap result left out
%        curr_stdev_one_bootstrap_left_out = std(curr_one_bootstrap_left_out_mean');  % stdev over results when one bootstrap sample is left out
% 
%        current_results = [current_results; curr_total_mean; curr_stdev_one_bootstrap_left_out; NaN .* ones(size(curr_total_mean))];
%        
% 
%    end
% 
%  end
% 
%  
% 
%  
% fprintf('\n') 
% disp(current_results)
% fprintf('\n')
% 


   
  



