function [correct_class_decision_values normalized_rank_results rank_confusion_matrix] = get_rank_and_decision_value_results(cv, YTe, classifier_labels, decision_values, create_rank_confusion_matrix)
% This helper method calculates the normalized rank results, the decision values for the correct class, and the rank confusion matrix information.
%  This method should be run every time the classifier is tested.
%
%  Input parameters:
%    YTe:  the labels for the current test points
%    decision_values:  The decision values for the current test point
%    create_rank_confusion_matrix:  whether a confusion matrix should be created
%
%  The returned values are in the structure NORMALIZED_RANK_AND_DECISION_VALUE_RESULTS and have the fields
%
%   correct_class_decision_values: a [num_test_points x 1] vector containing the decision value for the ith test point (for the correct class of the ith point)
%
%   normalized_rank_results:  a [num_test_points x 1] vector containing the normalized rank value for the ith test point 
%
%   rank_confusion_matrix: [num_predicted_classes x num_actual classes] current rank confusion matrix.  
%      the i, j entry of this matrix tells how high up on the rank of predictions for the jth class
%      was entry i.  
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


rank_confusion_matrix = [];


YTe_unique_values = unique(YTe);

    
% % %[vals sorted_inds] = sort(decision_values', 'descend');   % NEED TO BE CAREFUL THERE IS NOT TIE IN THE MAX DECISION VALUE (could create a bias in results)                        
[vals sorted_inds] = sort((decision_values + eps .* rand(size(decision_values)))', 'descend');   % adding a small amount of noise so that there is not a tie in the decision values                       
the_ranks_all_test_points = sorted_inds';   % each row has all the ranked order results for a given test point (i.e., all ranks for YTe(i))


% go through all test points and find the ranking and the decision value of the real label 
for iTestPoint = 1:length(YTe)                                            
    curr_rank_results(iTestPoint) = find(YTe_unique_values(the_ranks_all_test_points(iTestPoint, :)) == YTe(iTestPoint));  % find rank of the real label YTe 
    correct_class_decision_values(iTestPoint) = decision_values(iTestPoint, (classifier_labels == YTe(iTestPoint)));   % find the decision value for real label
end    



normalized_rank_results = 1 - (((curr_rank_results) - 1)./(length(YTe_unique_values) - 1));



% get information to create the rank confusion matrix
if create_rank_confusion_matrix == 1
    for iUniqueYTe = 1:length(YTe_unique_values)   % assuming that each time unique(YTe) has all possible label values (if it doesn't then this code will not work)

            curr_rank_inds = find(YTe == YTe_unique_values(iUniqueYTe));
            [rank_cm_vals rank_cm_inds] = sort(the_ranks_all_test_points(curr_rank_inds, :)');   % rank_cm_inds contains predicted classes from the most likely (first index value), to least likely (last index value) 

            rank_confusion_matrix(:, iUniqueYTe) = sum(1 -((rank_cm_inds -1)./(length(YTe_unique_values) - 1)), 2);  % add all test points together for a total rank confusion matrix

    end  
end


    
    





