function [all_XTr_pseudo all_XTe_pseudo] = turn_all_simultaneous_data_into_pseudo_populations(all_XTr_simul, all_YTr, all_XTe_simul, all_YTe, feature_to_siteID_mapping)

%  A helper function that takes simultaneously recorded data and turns it into pseudo-populations.  
%  This function should return results that are very similar to using normal pseudo-populations
%  (although this code will run much slower), but I include it anyway as a sanity check.  
%  The difference between this code and creating regular pseduo-populations is that 
%  when there are more trials recorded than are used in cross-validation splits, 
%  whole simultaneous trials from all neurons will be left out of the training and 
%  test data for this code, while regular pseudo-populations will take data from
%  all trials.
%
%  Notes about the code:
%   1.  Care is taken so that if population vectors contain neuron's responses from multiple time periods, 
%        when the data is shuffled, all data from the same trial (from a given site) will be permuted together.
%        This is important because mixing responses from the same trial (and site) into different population vectors could 
%        lead to biased results since there are dependencies between time points in a trial.  
%
%   2.  This code is slow.  It might be possible to write a faster version of the code (particularly 
%        if one doesn't have totally different randomization for each label, etc).  
%
%  Inputs:
%    all_XTr_simul{iTimePeriod}{iCV} = [num_features x num_points]  The simultaneously population training data
%    all_YTr The training labels  (can be a cell array for each time period or a vector having the same labels for all time periods)
%    all_XTe_simul{iTimePeriod}{iCV} = [num_features x num_points]  The simultaneously population test data
%    all_YTe The test labels  (can be a cell array for each time period or a vector having the same labels for all time periods)
%    feature_to_siteID_mapping{iTimePeriod} = a vector (length of num_features) that lists which site an a given feature 
%      in all_XTr_simul came from (this is important so that all features from a given site are moved together).  
%
%  Output:  all_XTr_pseudo, which is the simultaneous training data converted into pseudo-populations  
%           all_XTe_pseudo, which is the simultaneous test data converted into pseudo-populations  

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



for iTimePeriod = 1:length(all_XTr_simul)
    for iCV = 1:length(all_XTr_simul{1})

        if iscell(all_YTr)
            all_curr_labels = [all_YTr{iTime}; all_YTe{iTime}];   
        else
            all_curr_labels = [all_YTr; all_YTe];   
        end
        
        unique_labels = unique(all_curr_labels);
        unique_site_numbers = unique(feature_to_siteID_mapping{iTimePeriod});

               
        all_XTr_pseudo{iTimePeriod}{iCV} = NaN .* ones(size(all_XTr_simul{iTimePeriod}{iCV}));
        all_XTe_pseudo{iTimePeriod}{iCV} = NaN .* ones(size(all_XTe_simul{iTimePeriod}{iCV}));
        
        
        all_curr_data = [all_XTr_simul{iTimePeriod}{iCV} all_XTe_simul{iTimePeriod}{iCV}];
        
        
        
        for iLabel = unique_labels'
            
            all_curr_label_inds = find(all_curr_labels == iLabel);

            if iscell(all_YTr)
                curr_label_inds_Tr = find(all_YTr{iTime} == iLabel);    % labels in the original data (so I know where to replace things)
                curr_label_inds_Te = find(all_YTe{iTime} == iLabel);
            else
                curr_label_inds_Tr = find(all_YTr == iLabel);    % labels in the original data (so I know where to replace things)
                curr_label_inds_Te = find(all_YTe == iLabel);
            end
                
            if iscell(all_YTr)    
                num_Tr_labels = length(find(all_YTr{iTime} == iLabel));
            else
                num_Tr_labels = length(find(all_YTr == iLabel));
            end
            
                                    
            % moving around chunks of data from each neuron so that data from all times in a given trial is moved together
            for iSite = unique_site_numbers
                           
                curr_rand_label_inds = all_curr_label_inds(randperm(length(all_curr_label_inds))); 
                curr_rand_label_inds_Tr = curr_rand_label_inds(1:num_Tr_labels);    % permutated labels from combined training and test data (so that create pseduo-populations from both sets)
                curr_rand_label_inds_Te = curr_rand_label_inds((num_Tr_labels + 1):end);
                
                site_feature_inds  = find(feature_to_siteID_mapping{iTimePeriod} == iSite);   % ok, can save space by only storing the first column of all_XTr_feature_neuron_i                          
                all_XTr_pseudo{iTimePeriod}{iCV}(site_feature_inds, curr_label_inds_Tr) = all_curr_data(site_feature_inds, curr_rand_label_inds_Tr);   %[num_neuron_features x num_trials]
                all_XTe_pseudo{iTimePeriod}{iCV}(site_feature_inds, curr_label_inds_Te) = all_curr_data(site_feature_inds, curr_rand_label_inds_Te);   %[num_neuron_features x num_trials]
            
            end 
        end
               

    end
end














