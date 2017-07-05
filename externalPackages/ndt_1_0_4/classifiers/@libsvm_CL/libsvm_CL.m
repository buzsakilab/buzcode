classdef libsvm_CL

%  libsvm_CL is a Support Vector Machine (SVM) classifier object (CL) 
%   that uses the LIBSVM package (i.e., this function is Neural Decoding Toolbox
%   wrapper for the LIBSVM package).  A support vector machine is a classifier
%   that learns a function f that minimize the hinge loss between the training
%   data and labels, while also applying a penalty to more complex f (the penalty 
%   is based on the norm of f in a RKHS).  The SVM has a parameter C, that controls
%   the trade off between a smaller empirical loss (i.e., smaller prediction error 
%   on the training set), a more complex function f.  Support vector machine can use
%   different kernels to create nonlinear decision boundaries (several are supported here).
%
%  SVMs are designed to work on binary classification problems.  In order to 
%   support multi-class classification, we use two different methods.  The first method is 
%   called 'all-pairs' and works by training separate classifiers for all pairs of labels
%   (i.e., if there are 100 different classes then nchoosek(100, 2) = 4950 different classifiers
%   are trained.  Testing the classifier in all-pairs involves having all classifiers
%   classify the test point, and then the class label is given to the class the was chosen most
%   often by the binary classifiers (in the case of a tie in the number of classes that won a contest
%   the class label is randomly chosen).  The decision values for all-pairs are the number 
%   of contests won by each class (for each test point).
%
%  The second multi-class method is called 'one-vs-all' classification
%   and works by training one classifier for each class (thus if there are 100 classes
%   there will be 100 classifiers) using data from one class as the positive labels,
%   and data from all the other classes as the negative labels.  A test point is 
%   then run through these different classifiers and the class that has the largest 
%   SVM prediction value is returned as the predicted label 
%   (i.e., SVMs create a function f(x) = y, and the class label is 
%   usually given as sign(y), however here we are comparing the actual y values to determine
%   the label).  The decision values here are the f(x) = y values returned by the SVM. Our limited 
%   tests have found all pairs is faster and gives slightly more accuracy results so it is the
%   default (although the decision values might be considered more crude).
%
%
%  Function options and default values: 
%
%  The following values can be set to change this classifiers behavior (here we 
%  are assuming that svm = libsvm_CL)
%
%  - svm.C: scalar  (default svm.C = 1). This is the 1/regularization constant that determines
%       the trade off between the fit to the training data and the amount of regularization/simplicity
%       of the learned function.  A larger value of C means more emphasis on a better fit to the training data.
%
%  - svm.kernel:  'linear',  'polynomial'/'poly', or 'gaussian'/'rbf' (default 'linear).  The type of kernel used
%       which controls the class of functions that classifier is built from.  
%
%       -  If svm.kernel = 'polynomial', then the following options must be set:
%           - svm.poly_degree: scalar (there is no default value so this must be set).  The degree of the of the polynomial.
%           - svm.poly_offset: scalar (default svm.poly_offset = 0);  a constant offset if a polynomntal kernel is used 
%
%       -  If svm.kernel = 'gaussian', then the following options must be set:
%           - svm.gaussian_gamma (there is no default value so this must be set).  This controls the fall off
%             of the radial basis function kernel (larger values mean a slower fall off).
%
%  - svm.additional_libsvm_options (default svm.additional_libsvm_options = '').  This allows one to add additional
%     control of the classifiers behavior using a string that is in LIBSVM format.  For more details see: 
%     http://www.csie.ntu.edu.tw/~cjlin/libsvm/     
%
%  - svm.multiclass_classificaion_scheme (default value is 'all_pairs').  This field can
%    be set to 'all_pairs' or 'one_vs_all' and determines if a one-vs-all or an all-pairs
%    multi-class classification scheme is used, as described above (for binary problems
%    both all-pairs and one-vs-all will return the same predicted labels, although the 
%    decision values will differ).  
%
% Like all CL objects, there are two main methods, which are:
%
%  1.  cl = train(cl, XTr, YTr) 
%       This method takes the training data (XTr, YTr) to train the support vector machine.
%        XTr and XTe are in the form [num_features x num_examples]
%        YTr is in the form [num_examples x 1]
% 
%  2.  [predicted_labels decision_values] = test(cl, XTe)
%        This method takes the test data and produces a [num_test_points x 1] vector of 
%         predicted labels based on the  function learned from the training set.  It also 
%         returns a matrix of dimension [num_test_points x num_classes] of decision values.  
%
%
%  Notes:  
%   1.  If there is a tie among the decision values, then one of tied the classes
%          is chosen randomly as the predicted label.
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
        C = 1;  % inverse of the regularization constant - higher values cause more emphasis to be place on the
                %  emperical loss (and cause the classifier use more complex functions, which could lead to overfitting).
        kernel = 'linear';  % the type of kernel used in the SVM.  Choices are 'linear', 'polynomial'/'poly', or 'gaussian'/'rbf'.
        poly_degree = [];  % if a polynomial classifier is used, the degree of the polynomial must be set
        poly_offset = 1;  % a constant offset if a polynomntal kernel is used (default is to use an offset of 1)
        gaussian_gamma = [];  % this controlls the falloff of the Gaussian kernel (larger values mean slower falloff = wider Gaussian functions)
        additional_libsvm_options = '';  % can set this to use additional libsvm options that are not explicitly supported by this wrapper
        multiclass_classificaion_scheme = 'all_pairs';   % options are 'one_vs_all' or 'all_pairs'.
  end

  
  properties (GetAccess = 'public', SetAccess = 'private')
        labels = [];  % all the unique labels for each class 
        model = [];  % the model learned from the training data 
  end
    
    

    methods 

        % constructor 
        function cl = libsvm_CL
        end
        
  
        
        function cl = train(cl, XTr, YTr)  

            % added sanity check
            if size(YTr, 2) ~= size(XTr, 2)  &&  size(YTr, 1) ~= size(XTr, 2) 
               error('Number of columns in YTr, and XTr must be the same (i.e., there must be one and exactly one label for each data point)') 
            end

            cl.labels = unique(YTr);

            
            % create the libsvm arguments for the kernel
            cl.kernel = lower(cl.kernel);
            if strcmp(cl.kernel, 'poly'), cl.kernel = 'polynomial'; end
            if strcmp(cl.kernel, 'rbf'), cl.kernel = 'gaussian'; end
            
            switch cl.kernel
              case 'linear'
                kernelstring = '-t 0 ';
              case 'polynomial'
                kernelstring = ['-t 1 -g 1 -r ' cl.poly_offset ' -d ' cl.poly_degree ' '];
              case 'gaussian'
                kernelstring = ['-t 2 -g ' cl.gaussian_gamma ' '];
            end
               

            cstring = ['-c ' num2str(cl.C) ' '];  %  create the libsvm arguments for the (inverse of the) regularization constant       
            basicstring = ['-s 0 '];   % use C-SVC  (which is the default anyway)            
            probstring = '';  %['-b 1 '];
            
            
            paramstring = [basicstring cstring kernelstring probstring cl.additional_libsvm_options];
            
            
            
            if strcmp(cl.multiclass_classificaion_scheme, 'all_pairs')
               
                cl.model = svmtrain2(YTr, double(XTr'), paramstring);  % renamed the libsvm function svmtrain2 because the Matlab bioinformatics toolbox now has a funciton named svmtrain
                %cl.model = svmtrain(YTr, double(XTr'), paramstring);  

                
            elseif strcmp(cl.multiclass_classificaion_scheme, 'one_vs_all')
               
                % do one vs. all training
               for iClass = 1:length(cl.labels)
                   
                   reformatted_test_data = [XTr(:, (YTr == cl.labels(iClass))) XTr(:, (YTr ~= cl.labels(iClass)))];
                   curr_labels = zeros(size(YTr));
                   curr_labels(1:length(find(YTr == cl.labels(iClass)))) = 1;
                   the_weight = length(find(YTr ~= cl.labels(iClass)))./length(find(YTr == cl.labels(iClass)));
                   cl.model{iClass} = svmtrain2(curr_labels, double(reformatted_test_data'), paramstring); 
                   %cl.model{iClass} = svmtrain(curr_labels, double(reformatted_test_data'), paramstring); 
                       
               end
            end
            

        end
            
        
        
 
        function [predicted_labels decision_values] = test(cl, XTe)
        
            junk_fake_xte_labels = ones(size(XTe,2),1);
            
            
            if strcmp(cl.multiclass_classificaion_scheme, 'all_pairs')
            
                [predic_label, junk_fake_accuracy, libsvm_format_all_pairs_decision_vals] = svmpredict(junk_fake_xte_labels, double(XTe'), cl.model);   % edited by Ethan to work with libsvm 2.86 (compared to 2.8)

                predicted_labels = predic_label(:,1);
                if (size(predic_label,2)>1)
                    libsvm_format_all_pairs_decision_vals = predic_label(:,2);
                end


                % Creating the decision values that are the number of 'all pairs' classifiers that are in favor of a given class 
                %  (rather than the f(x) values that classifier returns)
               
                % For each test point, all of the 'all pairs' decision values are in each row of the libsvm_format_all_pairs_decision_vals matrix
                % the order goes:  [(1 vs 2), (1 vs 3), (1 vs 4), ... (1 vs num_classes), (2 vs 3), (2 vs 4)... (2 vs num_classes), ... etc... (num_classes -1 vs num_classes)]
                % Below I reparse this into a matrix of the form:  allpairs_decision_vals(iClass1, iClass2, iTestPoint) = libsvm_decision_value for iTestPoint destinguishing between iClass1 and iClass2
                % more negative values indicate preference for iClass1 and positive values indicate iClass2

                num_labels = length(cl.model.Label);
                end_inds = cumsum((cl.model.nr_class -1):-1:1);  
                start_inds = end_inds + 1;
                start_inds = [1 start_inds(1:end-1)];

                start_inds(end + 1) = 2;
                end_inds(end + 1) = 1;


                allpairs_decision_vals(1, :, :) = libsvm_format_all_pairs_decision_vals(:, start_inds(1):end_inds(1))';
                for iLabel = 2:num_labels

                    for iPrevLabel = 1:(iLabel -1)
                        curr_decision_vals(iPrevLabel, :) = squeeze(allpairs_decision_vals(iPrevLabel, iLabel-1, :) .* -1);
                    end
                    allpairs_decision_vals(iLabel, :, :) = [curr_decision_vals; libsvm_format_all_pairs_decision_vals(:, start_inds(iLabel):end_inds(iLabel))'];

                end

                % the final decision values are the number of all-pairs wins that a given class has
                for iPoints = 1:size(allpairs_decision_vals, 3)
                    decision_values(iPoints, :) = (sum(sign((allpairs_decision_vals(:, :, iPoints))), 2)');
                end

                
                % LibSVM's final result is based on the class that has the most 'all pairs wins' (and for ties it takes the lowest class number)
                % since this can introduce a bias (particularly in the confusion matrix), I will randomly select between classes that are tied by
                % adding a small amount of noise to the number of all pairs wins (I am doing this instead of using the randmax function b/c I want
                % this noise in the decision_values when the rank results and other functions are calculated from the decision values)
                
                decision_values = decision_values + eps .* rand(size(decision_values));
                
                [vals ind] = max(decision_values');
                predicted_labels = ind';

                predicted_labels = cl.labels(predicted_labels);
                
                
                
   
            elseif strcmp(cl.multiclass_classificaion_scheme, 'one_vs_all')
               
                
               % do one vs. all testing
               for iClass = 1:length(cl.labels)
                   
                   [predic_label, junk_fake_accuracy, libsvm_format_all_pairs_decision_vals] = svmpredict(junk_fake_xte_labels, double(XTe'), cl.model{iClass});
                   decision_values(:, iClass) = libsvm_format_all_pairs_decision_vals; 
            
                   % % linear f(x)  % the decision values are the same as these when a linear kernel is used...
                   %w = cl.model{iClass}.SVs' * cl.model{iClass}.sv_coef;
                   %linear_pred_vals(:, iClass) = (XTe' * w) - cl.model{iClass}.rho;
                   %all_weights(:, iClass) = w;
                   %all_b(:, iClass) = -1 .* cl.model{iClass}.rho;
               end
               
               
               [val ind] = randmax(decision_values');  
               predicted_labels = ind';
               
               predicted_labels = cl.labels(predicted_labels);
               
                              
            end   % end one vs all
            
            
        end  % end test method
        
    end  % end public methods
       
   
    
    

end   % end classdef








