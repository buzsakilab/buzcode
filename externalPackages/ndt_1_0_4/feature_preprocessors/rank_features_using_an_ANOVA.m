function [sorted_features ANOVA_pvals] = rank_features_using_an_ANOVA(X_data, Y_labels)

% A helper function that does a one-way ANOVA on each feature and returns a p-value
% for each feature and a ranking of the features based on their p-values.  
% Essentially this is a the same as running anova1 separately on each feature
% although this code is much faster.  
%
% Input:  
%       X_data is in the form [num_features x num_examples]
%       Y_labels is [num_examples x 1]   
%
% Output: 
%   sorted_features:  the features ranked from the smallest p-value (most significant feature) to the largest p-value (least significant feature)
%   ANOVA_pvals:  p-values for each feature  (in the original order of the input features).
%
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




unique_labels = unique(Y_labels);

for iLabel = 1:length(unique_labels)
    num_unique_labels(iLabel) = length(find(Y_labels == unique_labels(iLabel)));
end


if sum(abs(diff(num_unique_labels))) ~= 0   % if there are a different number of examples for some of the labels, we need to run a slightly slower version of the code
        
    for iLabel = 1:length(unique_labels)
        X_data_ANOVA_format{iLabel} = (X_data(:, (Y_labels == unique_labels(iLabel))))';   % use transpose so that we have num_points_in_group x num_features
    end

    ANOVA_pvals = oneway_ANOVA_unbalanced(X_data_ANOVA_format);

    
else     % all the classes have the same number of examples, so we can use a sligttly faster version of the code
    
    X_data_reshaped = zeros(num_unique_labels(1), length(unique_labels), size(X_data, 1));   % save time by preallocating memory
    
    for iLabel = 1:length(unique_labels)       
        X_data_reshaped(:, iLabel, :) = (X_data(:, (Y_labels == unique_labels(iLabel))))';
    end
    
    
    ANOVA_pvals = oneway_ANOVA_balanced(X_data_reshaped);

end


ANOVA_pvals = ANOVA_pvals';

[sorted_pvals, sorted_features] = sort(ANOVA_pvals);
    






function pvals = oneway_ANOVA_balanced(ANOVA_data)
% A helper function that calcuate the balanced one-way ANOVA p-value for each feature,
%   i.e., calculates the one-way ANOVA p-values when there are the SAME number of points 
%   in all groups. Since there is almost always the same number of examples for each 
%   class in the training data, I don't think the unbalanced version is actually needed
%   but I am keeping it here just in case.  
%   
% ANOVA_data is a tensor of:  [num_points_per_group x num_groups x num_features]


[num_points_per_group, num_groups, num_features] = size(ANOVA_data);


overall_mean = mean(mean(ANOVA_data, 1), 2);  % mean for all of the features (over all points)
mean_for_each_group = mean(ANOVA_data, 1);


between_group_degrees_of_freedom = (num_groups -1);
within_group_degrees_of_freedom = num_groups .* (num_points_per_group -1);    % (obvious change this if there are are different number of points per group)

      
RSS = num_points_per_group .* sum( (mean_for_each_group - repmat(overall_mean,[1,num_groups])).^2 ,2);   %  sum of how much the groups deviate from the grand mean (times of number of points per group)
TSS = sum(sum((ANOVA_data - repmat(overall_mean,[num_points_per_group, num_groups])).^2,1));  % total of how much all the points differ from the overall mean
SSE = TSS - RSS;    % difference between the total deviations from the group deviations


if (within_group_degrees_of_freedom > 0)
   MSE = SSE/within_group_degrees_of_freedom;
else
   MSE = NaN;
end

%F = repmat(Inf,[1,1,num_features]);     % should this be Inf's or zeros (want p-values of 1 for the case when mean for both classes is zero)
F = zeros(num_features, 1);
pvals = ones(num_features, 1);            

indices = find(SSE~=0);   % to prevent errors ...
F(indices) = (RSS(indices)/between_group_degrees_of_freedom) ./ MSE(indices);  

pvals(indices) = 1 - fcdf(F(indices), between_group_degrees_of_freedom, within_group_degrees_of_freedom);     % probability of the F ratio if the means of all groups are equal


pvals = squeeze(pvals);





function pvals = oneway_ANOVA_unbalanced(ANOVA_data)
% A helper function that calcuate the unbalanced ANOVA p-value for each feature, 
%   i.e., calculates a one-way ANOVA p-value when there are the DIFFERENT numbers of 
%   points in at least some of the groups.
%   ANOVA_data is a cell array:  ANOVA_data{num_groups} = [num_points_in_group x num_features]


num_groups = length(ANOVA_data);
num_features = size(ANOVA_data{1}, 2);   

for iGroup = 1:num_groups
    num_points_per_group(iGroup) = size(ANOVA_data{iGroup}, 1);    
    mean_for_each_group(iGroup, :) = mean(ANOVA_data{iGroup}, 1);
    
    sum_of_each_group(iGroup, :) = sum(ANOVA_data{iGroup}, 1);
    
    
    sum_of_squared_points_for_each_group(iGroup, :) = sum(ANOVA_data{iGroup}.^2, 1);   % used for alterative way to calculate p-values
    
end


overall_mean = sum(sum_of_each_group, 1)./repmat(sum(num_points_per_group), [1, num_features]);     % a real average over all points (groups with more points have a larger weight)


between_group_degrees_of_freedom = (num_groups -1);
within_group_degrees_of_freedom = sum(num_points_per_group - 1);    


for iGroup = 1:num_groups
    RSS_by_group(iGroup, :) = num_points_per_group(iGroup) .* (mean_for_each_group(iGroup, :) - overall_mean).^2;   %  sum of how much the groups deviate from the grand mean (times of number of points per group)
    TSS_by_group(iGroup, :) = sum((ANOVA_data{iGroup} - repmat(squeeze(overall_mean'), [1 num_points_per_group(iGroup)])').^2,1);  % total of how much all the points differ from the overall mean
end

RSS = sum(RSS_by_group, 1);  % between group sum of squares (SSb)  (also called treatment sum of squares)
TSS = sum(TSS_by_group, 1);  % total sum of squares (SStot)
SSE = TSS - RSS;    % difference between the total deviations from the group deviations (is equivalent to SSw since SStot = SSb + SSw)


if (within_group_degrees_of_freedom > 0)
   MSE = SSE/within_group_degrees_of_freedom;
else
   MSE = NaN;
end


F = zeros(num_features, 1);
pvals = ones(num_features, 1);            

indices = find(SSE~=0);   % to prevent errors ...
F(indices) = (RSS(indices)/between_group_degrees_of_freedom) ./ MSE(indices);  

pvals(indices) = 1 - fcdf(F(indices), between_group_degrees_of_freedom, within_group_degrees_of_freedom);     % probability of the F ratio if the means of all groups are equal


pvals = squeeze(pvals);





