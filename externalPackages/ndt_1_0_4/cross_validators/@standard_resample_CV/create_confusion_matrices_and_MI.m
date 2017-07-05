function DECODING_RESULTS = create_confusion_matrices_and_MI(cv, YTe, all_predicted_labels, DECODING_RESULTS)
%  A method to create confusion matrices for the zero-one loss and normalized rank results.
%   A confusion matrix is a matrix in which the columns express the true class label
%   and the rows express the predicted label.  For example, if column 2, and row 3
%   had a value of 7, it would mean that there were 7 times in which class 2 was mistakenly
%   predicted to be class 3 (summing the columns will give the total number of examples for 
%   each class, which can be used to turn the confusion matrix into a misclassification probability distribution).  

% confusion matrix => predicted_label x real_label x num_time_periods x (num_time_periods)
%
% I am assuming YTe is the same for all time periods, resampling iterations and CV splits (which should always be met as far as I can tell)


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



    real_labels = YTe;
    unique_labels = unique(YTe);

    % remap labels if they are not consequative positive numbers starting at 1
    if min(YTe) < 0 || ~(unique_labels(end) == length(unique_labels))

        old_all_predicted_labels = all_predicted_labels;
        all_predicted_labels = zeros(size(old_all_predicted_labels));

        for iNewLabels = 1:length(unique_labels)
            all_predicted_labels(find(old_all_predicted_labels == unique_labels(iNewLabels))) = iNewLabels;
            real_labels(find(YTe == unique_labels(iNewLabels))) = iNewLabels;
            cf_label_remapping(iNewLabels) = unique_labels(iNewLabels);  
        end            

        if length(find(all_predicted_labels == 0)) > 1
            warning('problems remapping labels, when creating the confusion matrix!!!!')
        end

         % save information about what original labels map on to the rows/columns of the confusion matrix
         DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.confusion_matrix_results.label_remapping = cf_label_remapping;  

    end   % end for remapping labels

    
    
    num_labels = length(unique_labels);
    num_train_time_intervals = size(all_predicted_labels, 4);   
    num_test_time_intervals = size(all_predicted_labels, 5);   
    num_predicted_labels = size(all_predicted_labels, 3);     
    num_resamples = size(all_predicted_labels, 1);
    
    

    
    % initialize variables the final variables that will store the confusion matrix and the separate test point confusion matrix
    %  (need to initialize these to zeros for the case when num_train_time_intervals ~= num_test_time_intervals) 
    if cv.confusion_matrix_params.create_confusion_matrix == 1
        if cv.confusion_matrix_params.save_confusion_matrix_only_train_and_test_at_same_time == 1
            confusion_matrix = zeros(num_labels, num_labels, num_train_time_intervals);
        else
            confusion_matrix = zeros(num_labels, num_labels, num_train_time_intervals, num_test_time_intervals);
        end
    end
    
    if cv.confusion_matrix_params.create_all_test_points_separate_confusion_matrix == 1
                
        if cv.confusion_matrix_params.save_confusion_matrix_only_train_and_test_at_same_time == 1
            all_test_points_separate_confusion_matrix = zeros(num_labels, num_predicted_labels, num_train_time_intervals);
        else
            all_test_points_separate_confusion_matrix =  zeros(num_labels, num_predicted_labels, num_train_time_intervals, num_test_time_intervals);      
        end
    end
    
    
    

    
   % begin the main loop that will create the confusion matrices and mutual information measures 
   for iTrain = 1:num_train_time_intervals
        for iTest = 1:num_test_time_intervals
            
            
            % initalize this variable to zero for each training and test time
            curr_CM = zeros(num_labels, num_labels);            
            if cv.confusion_matrix_params.create_all_test_points_separate_confusion_matrix == 1
                curr_CM_separate_for_each_test_point =  zeros(num_labels, num_predicted_labels);
            end
    
            
            
            for iResample = 1:num_resamples
                    
                
                % initalize this variable to zero each time the resample loop is run
                curr_resample_CM = zeros(num_labels, num_labels);

                                      
                for iCV = 1:size(all_predicted_labels, 2)

                    % put entry into into the confusion matrix:  
                    for iPredlabels = 1:num_predicted_labels

                        
                         % a [num_labels x num_labels] CM for the current resample run, train time and test time
                        curr_resample_CM(all_predicted_labels(iResample, iCV, iPredlabels, iTrain, iTest), real_labels(iPredlabels)) ...
                                = curr_resample_CM(all_predicted_labels(iResample, iCV, iPredlabels, iTrain, iTest), real_labels(iPredlabels)) + 1; 
                                                
                       % a [num_test_points x num_test_points] CM summed over all resample runs, for the current train time and test time                            
                        if cv.confusion_matrix_params.create_all_test_points_separate_confusion_matrix == 1
                                curr_CM_separate_for_each_test_point(all_predicted_labels(iResample, iCV, iPredlabels, iTrain, iTest), iPredlabels)...
                                    = curr_CM_separate_for_each_test_point(all_predicted_labels(iResample, iCV, iPredlabels, iTrain, iTest), iPredlabels) + 1;   
                        end
                                
                    end
                end


                
                % calculate the mutual information separately for each resample 
                if cv.save_results.mutual_information == 1  
                    prob_predicted_actual = curr_resample_CM./sum(sum(curr_resample_CM));
                    prob_predicted = sum(curr_resample_CM')./sum(sum(curr_resample_CM));
                    prob_actual = sum(curr_resample_CM)./sum(sum(curr_resample_CM));

                    MI_pieces_each_pred_actual = prob_predicted_actual .* log2(prob_predicted_actual .* 1./(prob_predicted' * prob_actual));
                    MI_pieces_each_pred_actual(find(isnan(MI_pieces_each_pred_actual))) = 0;  
                    DECODING_RESULTS.MUTUAL_INFORMATION.from_separate_confusion_matrix_for_each_resample.decoding_results(iResample, iTrain, iTest) = sum(sum(MI_pieces_each_pred_actual));
                end
                
                
                curr_CM = curr_CM + curr_resample_CM;    % curr_CM is the confusion matrix summed over all resample runs (and CV splits) - for the current train and test time
                
                clear curr_resample_CM   % actually need to set these to zeros rather than clearing them (but just clear them as a double sanity check)

                
            end   % end for resample runs
                
            
                       
            
            % if actually saving the full the confusion matrix
            if cv.confusion_matrix_params.create_confusion_matrix == 1 || cv.confusion_matrix_params.create_all_test_points_separate_confusion_matrix == 1

                % if only saving confusion matrix at the same training and test times and there are the same number of training and test times   
                if cv.confusion_matrix_params.save_confusion_matrix_only_train_and_test_at_same_time == 1 && (num_train_time_intervals == num_test_time_intervals) 

                    if iTrain == iTest   
                        confusion_matrix(:, :, iTrain) = curr_CM;
                        if cv.confusion_matrix_params.create_all_test_points_separate_confusion_matrix == 1
                            all_test_points_separate_confusion_matrix(:, :, iTrain) = curr_CM_separate_for_each_test_point;
                        end
                    end    


                % in the case when there are different number of training and test times, but save_confusion_matrix_only_train_and_test_at_same_time == 1 
                %  have each training time have the results summed over all test times available  
                elseif cv.confusion_matrix_params.save_confusion_matrix_only_train_and_test_at_same_time == 1  && (num_train_time_intervals ~= num_test_time_intervals) 

                    confusion_matrix(:, :, iTrain) = confusion_matrix(:, :, iTrain) + curr_CM;        
                    if cv.confusion_matrix_params.create_all_test_points_separate_confusion_matrix == 1
                      all_test_points_separate_confusion_matrix(:, :, iTrain) = all_test_points_separate_confusion_matrix(:, :, iTrain) + curr_CM_separate_for_each_test_point;
                    end

                    
                else

                    confusion_matrix(:, :, iTrain, iTest) = curr_CM;
                    if cv.confusion_matrix_params.create_all_test_points_separate_confusion_matrix == 1
                        all_test_points_separate_confusion_matrix(:, :, iTrain, iTest) = curr_CM_separate_for_each_test_point;
                    end

                end
            
                
            end


            
            % creating mutual information combined over all resample runs 
            if cv.save_results.mutual_information == 1  

                prob_predicted_actual = curr_CM./sum(sum(curr_CM));
                prob_predicted = sum(curr_CM')./sum(sum(curr_CM));
                prob_actual = sum(curr_CM)./sum(sum(curr_CM));

                MI_pieces_each_pred_actual = prob_predicted_actual .* log2(prob_predicted_actual .* 1./(prob_predicted' * prob_actual));
                MI_pieces_each_pred_actual(find(isnan(MI_pieces_each_pred_actual))) = 0;            
                DECODING_RESULTS.MUTUAL_INFORMATION.from_combined_confusion_matrix_over_all_resamples.decoding_results(iTrain, iTest) = sum(sum(MI_pieces_each_pred_actual));

            end

            
            clear curr_CM curr_CM_separate_for_each_test_point  % actually need to set these to zeros rather than clearing them (but just clear them as a double sanity check)
            
            
        end    % for iTest
    
   end   % for iTrain

   
   
   
   
  % save all the results 
   
    if cv.confusion_matrix_params.create_confusion_matrix == 1   % could avoid using this temporary variable, but not going to bother
        DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.confusion_matrix_results.confusion_matrix = confusion_matrix;
    end

    if cv.confusion_matrix_params.create_all_test_points_separate_confusion_matrix == 1
        DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.confusion_matrix_results.all_test_points_separate_confusion_matrix = all_test_points_separate_confusion_matrix;   % might be a bit memory intensive though
    end
   
               
   if cv.save_results.mutual_information == 1  
        DECODING_RESULTS.MUTUAL_INFORMATION.from_separate_confusion_matrix_for_each_resample.mean_decoding_results = squeeze(mean(DECODING_RESULTS.MUTUAL_INFORMATION.from_separate_confusion_matrix_for_each_resample.decoding_results, 1))';
        DECODING_RESULTS.MUTUAL_INFORMATION.from_separate_confusion_matrix_for_each_resample.stdev.over_resamples = squeeze(std(DECODING_RESULTS.MUTUAL_INFORMATION.from_separate_confusion_matrix_for_each_resample.decoding_results, [], 1))';
   end

   
   
   
   
   
   
% % a slow way to compute mutual information    
%curr_total = 0;
%for iPred = 1:size(curr_CM, 1)
%    for iActual = 1:size(curr_CM, 2)
%        curr_total = curr_total + prob_predicted_actual(iPred, iActual) .* log2(   prob_predicted_actual(iPred, iActual)./ (prob_predicted(iPred) .*   prob_actual(iActual)));
%    end
%end
% DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mutual_information(iTrain, iTest) = curr_total;   


     














