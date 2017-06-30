function DECODING_RESULTS = save_more_decoding_measures(cv, DECODING_RESULTS)
% A private method for standard_resample_CV that calculates additional decoding results and stdevs values
%  based on the original decoding results.  This function is called at the end of the decoding 
%  procedure (since these results are derived from other results).  

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




    DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results = squeeze(mean(mean(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.decoding_results, 2), 1));
    DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples = squeeze(std(mean(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.decoding_results, 2), [], 1));   

    
    if cv.save_results.normalized_rank == 1
        DECODING_RESULTS.NORMALIZED_RANK_RESULTS.mean_decoding_results = squeeze(mean(mean(DECODING_RESULTS.NORMALIZED_RANK_RESULTS.decoding_results, 2), 1));
        DECODING_RESULTS.NORMALIZED_RANK_RESULTS.stdev.over_resamples = squeeze(std(mean(DECODING_RESULTS.NORMALIZED_RANK_RESULTS.decoding_results, 2), [], 1));
    end

    if cv.save_results.decision_values == 1
        DECODING_RESULTS.DECISION_VALUES.mean_decoding_results = squeeze(mean(mean(DECODING_RESULTS.DECISION_VALUES.decoding_results, 2), 1));
        DECODING_RESULTS.DECISION_VALUES.stdev.over_resamples = squeeze(std(mean(DECODING_RESULTS.DECISION_VALUES.decoding_results, 2), [], 1));
    end


    if cv.save_results.ROC_AUC == 1   

         DECODING_RESULTS.ROC_AUC_RESULTS.separate_CV_ROC_results.mean_decoding_results = squeeze(mean(mean(mean(DECODING_RESULTS.ROC_AUC_RESULTS.separate_CV_ROC_results.decoding_results, 3), 2), 1));  
         DECODING_RESULTS.ROC_AUC_RESULTS.separate_CV_ROC_results.stdev.over_resamples = squeeze(std(mean(mean(DECODING_RESULTS.ROC_AUC_RESULTS.separate_CV_ROC_results.decoding_results, 3), 2), [], 1)); 

         DECODING_RESULTS.ROC_AUC_RESULTS.combined_CV_ROC_results.mean_decoding_results = squeeze(mean(mean(DECODING_RESULTS.ROC_AUC_RESULTS.combined_CV_ROC_results.decoding_results, 2), 1)); 
         DECODING_RESULTS.ROC_AUC_RESULTS.combined_CV_ROC_results.stdev.over_resamples = squeeze(std(mean(DECODING_RESULTS.ROC_AUC_RESULTS.combined_CV_ROC_results.decoding_results, 2), [], 1));       

    end


    for iTrain = 1:size(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.decoding_results, 3)
        for iTest = 1:size(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.decoding_results, 4)    

            temp = squeeze(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.decoding_results(:, :, iTrain, iTest)); 
            DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_CVs_combined_over_resamples(iTrain, iTest) = std(temp(:));

            if cv.save_results.normalized_rank == 1
                temp = squeeze(DECODING_RESULTS.NORMALIZED_RANK_RESULTS.decoding_results(:, :, iTrain, iTest)); 
                DECODING_RESULTS.NORMALIZED_RANK_RESULTS.stdev.over_CVs_combined_over_resamples(iTrain, iTest) = std(temp(:));
            end

            if cv.save_results.decision_values == 1
                temp = squeeze(DECODING_RESULTS.DECISION_VALUES.decoding_results(:, :, iTrain, iTest)); 
                DECODING_RESULTS.DECISION_VALUES.stdev.over_CVs_combined_over_resamples(iTrain, iTest) = std(temp(:));
            end

        end
    end

    



