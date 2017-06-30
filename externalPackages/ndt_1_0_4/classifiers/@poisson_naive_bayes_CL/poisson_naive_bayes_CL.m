classdef poisson_naive_bayes_CL

% poisson_naive_bayes_CL is a classifier (CL) object that implements a 
%   Poisson Naive Bayes classifier.  For each class, the expected number of occurances 
%   (denoted lambda) is calculated separately for each feature/neuron, by taking the mean
%   values from the training data for each class.  To evaluate whether a given test point 
%   belongs to class i, the log of the likelihood function is calculated using the lambda values as
%   parameters of Poisson distributions for each feature/neuron (i.e., a separate Poisson distribution 
%   is calculated for each neuron), and the probability of observing a given number of spikes for
%   a particular neuron is calculated using the lambda value for that 
%   neuron.  The overall likelihood value is calculated by multiplying the 
%   probabilities for each neuron together (i.e., Naive Bayes classifiers assume that each
%   neuron is independent), or equivalently, adding the log of the probabilities for
%   each neuron together.  The class with the highest likelihood value is chosen as the 
%   predicted label, and the decision values are the log likelihood values.
%
%
% Like all CL objects, there are two main methods, which are:
%
%  1.  cl = train(cl, XTr, YTr) 
%       This method takes the training data (XTr, YTr) and learns a mean vector
%       (i.e., the lambda values) for each class.
% 
%  2. [predicted_labels decision_values] = test(cl, XTe)
%       This method takes the test data and calculates the log of the likelihood
%       function for each class (which are the decision values), and returns the class
%       with the highest likelihood as the predicted_label.  
%
%  Notes:    
%
%  1. This classifier assumes that all feature values are integers.
%
%  2. If the estimate rate parameter (lambda) for any feature is 0 (i.e., if the mean of a feature in the training data is 0), 
%       then this 0 lambda estimate will be replaced by the value 1/(num_training_points_in_class_k +1); i.e., we will assume
%       that there is one more training point in which an event occurred.  This prevents errors if an event occurred on a test point
%       but the estimated lambda was 0.
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
        lambdas = [];   % the expected number of occurances for each neuron for each class
        labels = [];  % the unique labels for the different classes
    end



    methods 

        % constructor 
        function cl = poisson_naive_bayes_CL
        end
        
        
        function cl = train(cl, XTr, YTr)  

            % added sanity check
            if size(YTr, 2) ~= size(XTr, 2)  &&  size(YTr, 1) ~= size(XTr, 2) 
               error('Number of columns in YTr, and XTr must be the same (i.e., there must be one and exactly one label for each data point)') 
            end

            % sanity check to make sure the training data only contains integers
           if sum(sum(abs(round(XTr) - XTr))) > 10.^-2  % (not having it be exactly zero to account for round off error)
               error(['The training data must only contain integers to use this classifier.  Make sure the data was loaded as integers (e.g., using the appropriate flag in ' ...
                   'basic_DS/generalization_DS datasources, and that only appropriate feature preprocessors were used (e.g., do not use the zscore_normalize_FP)'])
           end
            
                
            
            cl.labels = unique(YTr);

            for iLabel = 1:length(cl.labels )
                
                lambdas(:, iLabel) = mean(XTr(:, (YTr == cl.labels(iLabel))), 2);
                
                % If the rate parameter (lambda) is zero, there can be problems if the test data has a value greater than zero
                % (i.e., will get zero probability, which creates problems when taking the log).
                % We will deal with this problem by assuming that there is one more data point with a 1 (when all training data has zero occurances)
                zero_lambda_replacement_value = 1/(length(find(YTr == cl.labels(iLabel))) + 1);   % if found zeros for all trials, assume that there is one more trial where a 1 was found...
                lambdas((lambdas(:, iLabel) == 0), iLabel) = zero_lambda_replacement_value;
                
            end

            cl.lambdas = lambdas;
            
        end
            
        
        
        
        function [predicted_labels decision_values] = test(cl, XTe)
            
            
            % sanity check to make sure the test data only contains integers
           if sum(sum(abs(round(XTe) - XTe))) > 10.^-2   % (not having it be exactly zero to account for round off error)
               error(['The test data must only contain integers to use this classifier.  Make sure the data was loaded as integers (e.g., using the appropriate flag in ' ...
                   'basic_DS/generalization_DS datasources, and that only appropriate feature preprocessors were used (e.g., do not use the zscore_normalize_FP)'])
           end
            
            
            
            % compute simultaneously all p-values for all features and test points (easy b/c everything is independent)
            curr_lambdas = repmat(cl.lambdas, [1, 1, size(XTe, 2)]);
            XTe_repmat_for_all_classes = permute(repmat(XTe, [1, 1, size(cl.lambdas, 2)]), [1 3 2]);
            %p_vals = exp(-curr_lambdas + XTe_repmat_for_all_classes  .* log(curr_lambdas) - gammaln(XTe_repmat_for_all_classes  + 1));
            %log_likelihoods = squeeze(sum(log(p_vals)));
            
            log_likelihoods = squeeze(sum(-curr_lambdas + XTe_repmat_for_all_classes  .* log(curr_lambdas) - gammaln(XTe_repmat_for_all_classes  + 1), 1)); % possibly even faster...
            % % %log_likelihoods = squeeze(sum(-curr_lambdas + XTe_repmat_for_all_classes  .* log(curr_lambdas))); % since the factorial term does not depend on the class, code might be faster
                                                                                                              % if this term is removed.  The decision values returned
                                                                                                              % will no longer be exactly the log likelihood values but rather the log likelihood
                                                                                                              % values minus a constant.  Using this will make AUROC values worse, so don't do it!!!!
            
           % prior to version 1.0.4 the following line was used.  I added a the sum over the first dimention (i.e., sum(XXX, 1)) so that the code will run even when only a single site is used                                                                                                   
           %log_likelihoods = squeeze(sum(-curr_lambdas + XTe_repmat_for_all_classes  .* log(curr_lambdas) - gammaln(XTe_repmat_for_all_classes  + 1))); % possibly even faster...
                                                                                                  
                                                                                                              
            [vals inds] = randmax(log_likelihoods);
            predicted_labels  = cl.labels(inds);
            decision_values = log_likelihoods';
            

        end
        
                
    end  % end public methods
       
    

       

end   % end classdef








