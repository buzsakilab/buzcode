classdef max_correlation_coefficient_CL

% max_correlation_coefficient_CL is a classifier object (CL) 
%  that learns a mean population vector (template) for each class from the 
%  training data.  When the classifier is tested, the correlation coefficient
%  is calculated between a test point and each of the class templates, and the 
%  class with the largest correlation coefficient value is selected as the label
%  (and the decision values are the corrcoef values).  If there is a only one
%  features in XTr (and XTe), the decision is made based on the shortest squared
%  deviation between training and test value (and the negative of these distances are 
%  returned as decision values). If two or more classes are equally well correlated 
%  with a test point, then one of the tied classes is chosen randomly as the predicted label.  
%
%
% Like all CL objects, there are two main methods, which are:
%
%  1.  cl = train(cl, XTr, YTr) 
%         This method takes the training data (XTr, YTr) and learns a mean vector
%           (i.e., a template) for each class.
% 
%  2.  [predicted_labels decision_values] = test(cl, XTe)
%         This method takes the test data and calculates the correlation coefficient
%           value between each test point and each learned class template.  The 
%           predicted label for a test point is the class that had the largest corrcoef
%           value with the test point, and the decision values are the actually corrcoef values.
%
%
%  Notes:  
%    1.  If there is only one feature in XTr and XTe, then the prediction 
%         is made based on the negative of the squared difference between the test
%         feature and the training feature for each class (and these are the
%         decision values that are returned).
%
%    2.  If there is a tie among the decision values  (i.e., if the corrcoef
%          value is the same for 2 different classes), then one of tied the classes
%          is chosen randomly as the predicted label.
%  
%
%   XTr and XTe are in the form [num_features x num_examples]
%   YTr is in the form [num_examples x 1]
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




  properties 
        templates = [];   % average of the training vectors for each class
        labels = [];  % all the unique labels for each class (one for each template)  
  end



    methods 

        % constructor 
        function cl = max_correlation_coefficient_CL
        end
        
        
        function cl = train(cl, XTr, YTr)  

            % added sanity check
            if size(YTr, 2) ~= size(XTr, 2)  &&  size(YTr, 1) ~= size(XTr, 2) 
               error('Number of columns in YTr, and XTr must be the same (i.e., there must be one and exactly one label for each data point)') 
            end

            unique_labels = unique(YTr);

            for i = 1:length(unique_labels)
                template(:, i) = mean(XTr(:, (YTr == unique_labels(i))), 2);
            end

            cl.templates = template;
            cl.labels = unique_labels;
            
        end
            

        
        function [predicted_labels decision_values] = test(cl, XTe)
        
                  
            if size(XTe, 1) > 1   % if each data point is only 1 dimensional (could return an error, but instead just returning the template that had the closest value)
           
                % % use the corrcoef function to get the correlation coefficients
                % all_correlation_coefficients = corrcoef([cl.templates, XTe]);  % comment this out if using Octave and uncomment the line below
                % template_corrcoeffs = all_correlation_coefficients(size(cl.templates, 2) +1:end, 1:size(cl.templates, 2));
                
                % another way to compute the correlation coefficients - switched my code to this to make it compatible with Octave
                 mean_subtracted_templates = cl.templates - repmat(mean(cl.templates), [size(cl.templates, 1) 1]);
                 mean_subtracted_XTe = XTe - repmat(mean(XTe), [size(XTe, 1) 1]);
                 normalization_matrix = sqrt(diag(mean_subtracted_templates' * mean_subtracted_templates)) * sqrt(diag(mean_subtracted_XTe' * XTe))';
                 template_corrcoeffs = ((mean_subtracted_templates' * mean_subtracted_XTe)./normalization_matrix)';

            
            else   %  if there is only one feature, select the class with closest value to that feature

                % the squared difference between each class mean and each test point  (which are both scalars)
                template_corrcoeffs  = -1 .* (repmat(cl.templates, [size(XTe, 2), 1]) - repmat(XTe', [1 size(cl.templates, 2)])).^2;   
                
            end

            
            [val ind] = randmax(template_corrcoeffs');   % using randmax to deal with ties in max correlation value
            predicted_labels = cl.labels(ind);
            decision_values = template_corrcoeffs; 


            if (size(template_corrcoeffs, 1) .* size(template_corrcoeffs, 2)  ~= sum(sum(isfinite(template_corrcoeffs))))
               warning('this matrix contains some numbers that are not finite!!!')
            end
            
            
        end
        
    end  % end public methods
       
   

   
end   % end classdef








