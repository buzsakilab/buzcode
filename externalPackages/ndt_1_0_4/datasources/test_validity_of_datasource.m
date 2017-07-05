function test_validity_of_datasource(ds)

% A helper function that tests whether a given data source returns any of the same data
%  in the training and test data.  Obviously if some of the same data is in the training 
%  and test sets then there are major problems with the datasource!
%
% The input of this function should be a datasource that has been initialized with 
%  real data.  This function replaces the real data with unique numbers and tests
%  whether any of the same numbers are in the training and test sets.  The function
%  prints a message that says whether the datasource has problems or seems ok.


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




% replace the original data with unique numbers for data point

the_original_data = ds.the_data;


start_ind = 1;

for iSite = 1:length(the_original_data)

    end_ind = start_ind - 1 + prod(size(the_original_data{iSite}));
    
    the_new_data{iSite} = reshape(start_ind:end_ind, [size(the_original_data{iSite}, 1) size(the_original_data{iSite}, 2)]);

    start_ind = end_ind + 1;

end


clear the_original_data

ds.the_data = the_new_data;



% get the training and test data with the new data that has unique labels at each time point, and 
%   make sure that no unique label is in both the training and test sets

[XTr_all_time_cv YTr_all_cv XTe_all_time_cv YTe_all_cv] = get_data(ds);


num_probs_with_ds = 0;


for iTime = 1:size(XTr_all_time_cv, 2)   
    for iCV = 1:ds.num_cv_splits   
        
        if length(intersect(XTr_all_time_cv{iTime}{iCV}(:), XTe_all_time_cv{iTime}{iCV}(:)))
            num_probs_with_ds = num_probs_with_ds + 1;
        end
                
    end
end


if num_probs_with_ds > 1
    'This datasource has problems - some of the same data is in the training and test sets!'
else
    'This datasource seems ok'
end



