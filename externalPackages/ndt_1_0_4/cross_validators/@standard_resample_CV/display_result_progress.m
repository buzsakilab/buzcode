function display_result_progress(cv, results_to_display, display_progress)

 % Displays the current mean and stdev when a single resample result is left out for different result types
 %  to give an idea whether the algorithm has converged to a stable decoding accuracy
 %  (these values could be used for stopping criteria, rather than running a fixed number of resampling runs)

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
 
 


i = 1;
num_training_times = 0;
while num_training_times(1) == 0
    if i ~= 4
        num_training_times = size(results_to_display{i}, 3);
        num_test_times = size(results_to_display{i},4);
    else
        num_training_times = size(results_to_display{i}, 4);
        num_test_times = size(results_to_display{i}, 5);
    end
    i = i + 1;
end



 % default behavior if there is only 1 training point and display_progress.training_time_to_display_results == -1, is to just use this one training time for the results
 if (display_progress.training_time_to_display_results == -1)  && (num_training_times == 1)
     display_progress.training_time_to_display_results = 1;
 end


 % display_progress.training_time_to_display_results should be set if the num_training_times ~= num_test_times
 if (display_progress.training_time_to_display_results == -1) && (num_training_times~= num_test_times) && (num_test_times ~= 1)

     warning(['If the number of training times does not equal the number of test times, and one wants to display the resampled CV progress, ' ...
                'then display_progress.training_time_to_display_results should be set.  Using the default value of display_progress.training_time_to_display_results = 1, ' ...
                'to display the results here.'])

     display_progress.training_time_to_display_results = 1;

 end
 
 
 
 
current_results = [];
 

 for iResult = 1:length(results_to_display)
       
              
     iResample = size(results_to_display{1}, 1);

     if iResample == 0
         continue;
     end
     

    % if a particular training time that has been specified to display the results
    % note: if only saving results when training and testing at the same time period, then the displayed results will only be for one time point 
    if display_progress.training_time_to_display_results > 0
          

        if iResult ~= 4
           curr_bootstrap_progress_results = squeeze(mean(results_to_display{iResult}(:, :, display_progress.training_time_to_display_results, :), 2));  
        else
           curr_bootstrap_progress_results = squeeze(mean(mean(results_to_display{iResult}(:, :, :, display_progress.training_time_to_display_results, :), 3), 2));  
        end
        
        if iResample == 1
           curr_bootstrap_progress_results = curr_bootstrap_progress_results';
        end
        
    else
        
        % if only saving the results for training and testing at the same time
        if ((length(size(results_to_display{iResult})) == 3) && (iResult ~=4)) || ((length(size(results_to_display{iResult})) == 4)  && (iResult == 4))
            
            if iResult ~=4
                 curr_bootstrap_progress_results = squeeze(mean(results_to_display{iResult}, 2));
            else
                 curr_bootstrap_progress_results = squeeze(mean(mean(results_to_display{iResult}, 3), 2));
            end
            
            if iResample == 1
               curr_bootstrap_progress_results = curr_bootstrap_progress_results';
            end            
            
        else
        
             for iDisplayResultTrainingTime = 1:num_training_times 

                 if iResult ~= 4
                    curr_bootstrap_progress_results(:, iDisplayResultTrainingTime) = squeeze(mean(results_to_display{iResult}(:, :, iDisplayResultTrainingTime, iDisplayResultTrainingTime), 2));
                 else
                    curr_bootstrap_progress_results(:, iDisplayResultTrainingTime)  = squeeze(mean(mean(results_to_display{iResult}(:, :, :, iDisplayResultTrainingTime, iDisplayResultTrainingTime), 3), 2));
                 end
             end
         
        end
         
         
         
    end
        
        
     
   if iResample == 1
       current_results = [current_results; curr_bootstrap_progress_results; NaN .* ones(size(curr_bootstrap_progress_results))];   % can't show cofficients of variation if only 1 bootstrap has been run (so skip the following part)
   else

       curr_total_mean = mean(curr_bootstrap_progress_results, 1);
       curr_one_bootstrap_left_out_mean = ((repmat(curr_total_mean, [iResample 1]) - (curr_bootstrap_progress_results./iResample)) .* (iResample./(iResample -1)))';   % means over bootstraps with the ith bootstrap result left out
       curr_stdev_one_bootstrap_left_out = std(curr_one_bootstrap_left_out_mean');  % stdev over results when one bootstrap sample is left out

       current_results = [current_results; curr_total_mean; curr_stdev_one_bootstrap_left_out; NaN .* ones(size(curr_total_mean))];
       

   end

 end

 

 
fprintf('\n') 
disp(current_results)
fprintf('\n')



   
  



