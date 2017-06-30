function all_XTr_pseudo = turn_training_simultaneous_data_into_pseudo_populations(all_XTr_simul, all_YTr, feature_to_siteID_mapping)

% A helper function takes simultaneously recorded training data and turns them into pseudo-populations.  
%  This is useful for assessing whether training on shuffled data leads to similar performance (and decision boundaries)
%  as training with simultaneously recorded data (i.e., this allows one to compute I_diag as discussed by Averbeck, Latham and Pouget, 
%  in Neural corelations, population coding, and computation, Nature Neuroscience May 2006).  
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
%    feature_to_siteID_mapping{iTimePeriod} = a vector (length of num_features) that lists which site an a given feature 
%      in all_XTr_simul came from (this is important so that all features from a given site are moved together).  
%
%  Output:  all_XTr_pseudo, which is the simultaneous data converted into pseudo-populations  

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
            curr_labels = all_YTr{iTimePeriod};
        else
            curr_labels = all_YTr;
        end
        unique_labels = unique(curr_labels);
        
        unique_site_numbers = unique(feature_to_siteID_mapping{iTimePeriod});

        all_XTr_pseudo{iTimePeriod}{iCV} = NaN .* ones(size(all_XTr_simul{iTimePeriod}{iCV}));
        
        
        for iLabel = unique_labels'
            
            curr_label_inds = find(curr_labels == iLabel);

            for iSite = unique_site_numbers
 
                curr_rand_label_inds = curr_label_inds(randperm(length(curr_label_inds))); 
                
                site_feature_inds  = find(feature_to_siteID_mapping{iTimePeriod} == iSite);    
                          
                all_XTr_pseudo{iTimePeriod}{iCV}(site_feature_inds, curr_label_inds) = all_XTr_simul{iTimePeriod}{iCV}(site_feature_inds, curr_rand_label_inds);   %[num_neuron_features x num_trials]
            end 
        end
               

    end
end


 








