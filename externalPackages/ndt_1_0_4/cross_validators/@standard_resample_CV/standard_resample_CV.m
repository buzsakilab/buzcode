classdef standard_resample_CV

% standard_resample_CV is an object that implements the main functionality of
%   a cross-validator (CV) object.  Namely, it takes a classifier (CL), a datasource (DS)
%   and optionally feature preprocessor (FP) objects, and it run a cross-validation
%   decoding scheme by training and testing the classifier with data generated
%   from the datasource object (and possibly applying feature pre-processing first).
%   The cross-validation procedure is run multiple times
%   in which the data is resampled so that the datasource generates different training
%   and test splits of the data, and the final decoding performance is averaged over these
%   resampled runs.  
%
% The run_cv_decoding method for this object returns a "DECODING_RESULTS" structure
%   that contains zero-one loss results as well as potentially other types of results
%   (including normalized rank results, average decision values, and area under ROC curve results).
%   Confusion matrices and different measures of the variability of the decoding results
%   can also be returned.  The options that allow one to change the behavior of this
%   object and the results that this object returns are described in more detail below.
%
% The constructor for this object has the form:  
%   cv = basic_resample_CV(the_datasource, the_classifier, the_feature_preprocessors)
%   where the_datasource is a DS object, the_classifier is a CL object, and 
%   the_feature_preprocessors is an optional cell array that contains FP objects.  
%
%
% Properties that can be set are:
%
%   num_resample_runs (default = 50), which is the minimum number of resample runs
%      that the decoding procedure will be run.  
%
%   test_only_at_training_times (default = 0).  If this is set to zero than the classifier
%      will be trained and testing at all times, which will create a temporal-cross-training (TCT) 
%      matrix of results that have the decoding accuracies for training at time 1 and testing at time 2.  
%      This allows one to test if the neural code is stationary (see Meyers et al, 2008).
%      If this properties is set to 1, then the results will only be for training and testing
%      at the same time - which can cause the run_cv_decoding method to potentially run 
%      much faster, particularly if it classifier is slow at returning test results.
%        
%
%   stop_resample_runs_only_when_specfic_results_have_converged. Setting the different fields in this structure causes 
%       the resample runs to continue beyond those specified by num_resample_runs if the results have not converged to a 
%       stable estimate. stop_resample_runs_only_when_specfic_results_have_converged has the following fields 
%       that can be set which control the resample run stopping criteria for different types of results:
%           .zero_one_loss_results: controls whether the zero one loss results have converged
%           .normalized_rank_results: controls whether the normalized rank results have converged
%           .decision_values: controls whether the decision value results have converged
%           .combined_CV_ROC_results: controls whether the combined CV ROC results have converged
%           .separate_CV_ROC_results: controls whether the separate CV ROC results have converged
%      By default all these fields are set to empty meaning that no converge stoping criteria are set by default.
%      Setting any of these fields to a particular value causes the run_resample_cv method 
%      to keep running the resample run loop until the given mean result (over resample runs) changes by less than 
%      the specified value (over all training and test time periods) when any one resample run is left out.
%      This can be useful for speeding up the run-time when the decoding results have converged to a stable value. For example, one could 
%      set num_resample_runs to a lower number (say 10) and then setting  .zero_one_loss_results to a smallish value (say 1), which might cause 
%      fewer than the default value of 50 resample runs to be executed while getting results that are almost as accurate - i.e., there would 
%      at most be a change of 1 in the decoding accuracy (and any point in time) if the most sensitive resample run was left out. If any of 
%      these fields are set, there will still be a minimum number of resample runs executed that is specified by the num_resample_runs property, 
%      and then there will be additional resample runs to be executed until the desired convergence level is achieved. There is also an
%      additional field .stop_criteria_is_absolute_result_value (default = 1), which specifies whether the value set should be taken as an 
%      absolute change in the decoding accuracy - e.g., the actual zero-one decoding result values should change by less than 1 when any
%      resample run is left out. If this field is set to 0, then the values specified are given as a percentage change that should not occur
%      if any resample run is left out relative to the maximum of the mean decoding acccuracy - i.e., a value of 1 would 
%      mean that the results of leaving the ith resample run out, should not chance by more than 1% at that time point relative to the maximum 
%      decoding accuracy achieved (since the scale of a plot is determined relative to the maximum decoding accuracy this shows how much variance there is
%      in the results on a plot due to not using more resample runs). 
%      
%      
%     
% If a classifier is used that returns decision values along with zero-one loss results,
%  then the additional decoding accuracy results can be calculated and returned by setting
%  the following options:
%
%   save_results    
%       .normalized_rank (default = 1), returns normalized rank results
%   	.decision_values (default = 1), returns the decision values
%       .extended_decision_values (default is 0). If this property is set to 1 
%           then all the decision values for the correct class are saved. If 
%           this property is set to 2, then all decision values from all classes are
%           saved (doing this can potentially take up a lot of memory/disk space).
%       .ROC_AUC (default = 1), returns area under ROC curve results.  These results 
%           are calculated separately for each class k, with the test decision values from 
%           class k being the positive examples, and all other decision values being the negative
%           examples.  The results can be calculated separately on the test points on each CV split, 
%           (save_results.ROC_AUC = 3) or from combining all decision values over all cross-validation splits
%           (save_results.ROC_AUC = 2).  If save_results.ROC_AUC = 1, then both separate CV results, and combined
%           CV results will be saved.  
%       .mutual_information (default = 1).  returns estimates of mutual information
%           calculated from the confusion matrix.
%           
%
%  A number of parameters can be set to create confusion matrices for the results.
%   A confusion matrix is a matrix in which the columns express the true class label
%   and the rows express the predicted label.  For example, if column 2, and row 3
%   had a value of 7, it would mean that there were 7 times in which data from class 2 was mistakenly
%   predicted to be from class 3 (summing the columns will give the total number of examples for 
%   each class, which can be used to turn the confusion matrix into a misclassification probability
%   distribution).  The following parameters can be set to create confusion matrices
%
%   confunsion_matrix_params
%       .create_confusion_matrix (defatul is 1).  Create a basic zero-one loss confusion matrix
%       .save_confusion_matrix_only_train_and_test_at_same_time (default is 1). Saves the 
%           confusion matrix only for when training and testing was done at the same time.
%       .create_all_test_points_separate_confusion_matrix (default is 0).  Creates a confusion matrix 
%           where all test points are separate, i.e., if there are two test points from the same
%           class on a given CV split, these two test points will be given separate column
%           in the confusion matrix (this is useful when ID are remapped using the 
%           generalization_DS datasource).
%
%                    
%  The display_progress options allow one to display the decoding result 
%   (for different result types) as the decoding procedure as it is running.
%   The results displayed show the mean decoding results, as well as as a measure
%   of the variability of the results (which gives a sense of whether enough resample
%   iterations have been run so that the results have converged to a stable solution). The measure of
%   variability of the results is calculated by computing the mean decoding results 
%   if the i-th resample run was left out (i.e., if there was one less resample run). This
%   is done for each resample run, and the standard deviation is taken over these 
%   one-resample-left-out means.  Thus this measure gives a sense of how much the results
%   would vary if one less resample iteration was run (one can get a rough sense of whether
%   the results have converged if one compares these to the mean results - overall these numbers
%   should be very small).  The following options allow one to display progress of different
%   result types (the results for different types are displayed for all times, and are separated
%   by NaNs, i.e,. mean_type1 stdev_type1 NaN mean_type2 stdev_type2, NaN, etc.).
%
%   display_progress
%       .resample_run_time (default = 1).  Will display how many resample runs have completed, the amount of
%           time the last resample run took, and an estimate of the time when the code will be done running.
%       .zero_one_loss (default = 1).  Display zero-one loss results         
%       .normalized_rank (default = 0). Display normalized rank results
%       .decision_values (default = 0). Display decision values
%       .separate_CV_ROC_results (default = 0).  Display ROC AUC results computed separately on each CV split
%       .combined_CV_ROC_results (default = 0). Display ROC AUC results combined points from all CV splits
%       .training_time_to_display_results (default = -1).  The training time for which to display the test result values. 
%          A value of -1 means that the displayed results will for training at the test time points 
%          (or if only 1 training time was used, then a value of -1 means use that time for all test times). 
%       .convergence_values (default = 0). Displays the current resample run convergence values.
%
%  The run_cv_decoding method returns a structure called DECODING_RESULTS that contains
%   the following properties.
%
% DECODING_RESULTS
%
%   .FP_INFO{iFP}   If any of the feature preprocessing algorithms had returned information to be
%       saved via their get_current_info_to_save method, then this structure will contain 
%       the additional information that the FP algorithm wanted to be saved.  
%
%   .CV_PARAMETERS  This structure contains additional information about parameters
%       used in the decoding procedure.  This structure has the following fields:
%
%       .num_test_points_on_each_CV_run      The number of test points on each CV split
%       .num_training_points_on_each_CV_run  The number of training points on each CV split
%       .dimension_of_data_points            The dimensionality of the training/test points
%       .unique_labels                       A vector containing the unique test labels that were used
%       .num_unique_labels                   The number of unique labels (classes) 
%       .num_CV_splits                       The number of cross-validation splits of the data
%       .num_resample_runs                   The number of resample runs used
%       .toolbox_version_number              The number of current neural decoding toolbox being used
%       .classifier_name                     The class name of the classifier used 
%       .feature_preprocessor_names          The class names of the feature preprocessors used 
%
%    .DS_PARAMETERS  If the datasource has a method called get_properties, then the properties
%       returned by this method will be saved in this structure.
%
%   .ZERO_ONE_LOSS_RESULTS  This structure contains the main (zero-one loss) decoding results.
%       The following results can be returned.
%
%       .decoding_results [num_resample_runs x num_CV_splits x num_training_times x num_test_times]
%           This tensor contains all the decoding results separately for each resample run and
%           cross-validation split.
%       .mean_decoding_results: [num_training_times x num_test_times]  This contains the mean
%           decoding results averaged over all resample runs and CV splits
%       .confusion_matrix_results  If the confusion matrix field options are set (e.g., 
%           cv.confusion_matrix_params.create_confusion_matrix == 1), then this structure
%           can contain the following fields:
%           .confusion_matrix [num_predicted_classes x num_actual_classes x num_training_times x num_test_times]
%               This contains the confusion matrix specifying how often a test point j was classified as 
%               belonging to class i (where j indicates column indecies, and i indicates row indecies).  If 
%               .save_confusion_matrix_only_train_and_test_at_same_time = 1, then the confusion matrix is
%               [num_predicted_classes x num_actual_classes x num_training_times] large and only contains the confusion
%               matrices when training and testing at the same time period (which saves a lot of disk space when saving the results).
%           .label_remapping  If the labels given to the classifier are not consecutive integers, this 
%               vector indicates how the labels have been remapped on to the columns of the confusion matrix
%               (i.e., the first value in the vector indicates the class in the first row of the confusion matrix, etc.).
%           .all_test_points_separate_confusion_matrix If the create_all_test_points_separate_confusion_matrix flag is set to 1, 
%               then this variable contains a confusion matrix that is [num_test_points x num_actual_classes x num_training_times x num_test_times]
%               large, with each test point in a given CV split being given a separate row in the confusion matrix (this confusion matrix will 
%               only differ from the regular confusion matrix if there are multiple test points from the same class in
%               a given CV split).  This matrix can be useful if the labels have been remapped to different classes using
%               the generalization_DS datasource in order to see the confusion in the original unremapped labels
%       .stdev  This structure contains information about the variability of the results.  This structure has the following fields:
%           .all_single_CV_vals(num_resmaples, num_CV, num_train, num_test)
%               Each CV run produces a value for each test point (i.e., 0 and 1's, 
%               normalized ranks, or decision values).  This matrix contains the stdevs
%               over these values for each CV run.
%           .all_single_CV_vals_combined(num_resmaples, num_train, num_test)
%               This is the same as .all_single_CV_vals but all the values from the
%               CV runs are combined together first before
%               the stdev is taken (this is done separately for each resample run).  
%           .over_CVs(num_resmaples, num_train, num_test)
%               The mean decoding results in each CV run are computed (i.e., the mean
%               of the 0, 1's, ranks or decision values of the points in a CV split)
%               and the stdev over these mean CVs are calculated (this is done separately
%               for each resample run).
%           .over_CVs_combined_over_resamples(num_train, num_test)
%               This is the same as .over_CVs but all the values from all the resample 
%               runs are combined before the stdev is taken.
%           .over_resamples(num_train, num_test)
%               This calculates the mean over all the (mean) CV values for a resample run,
%               and then calculates the stdev over the different resample runs.  
%
%         It should be noted that .over_CVs, .over_CVs_combined_over_resamples and
%           .over_resamples could have all been computed running the decoding experiment
%           using the values in .decoding_results, but we precompute them here for convenience.
%
%
%    For classifiers that return decision values, this CV object can return additional results in 
%     the structures NORMALIZED_RANK_RESULTS, DECISION_VALUES and ROC_AUC_RESULTS that are described below.  
%
%    .NORMALIZED_RANK_RESULTS  This structure contains the normalized rank results.  
%       Normalized rank results are the results based on using the prediction values from
%       the classifier to create an ordered list of predictions (i.e., the most likely class is x, 
%       second most likely class is y, etc.), and then assessing how far down on the list is the
%       correct label.  The results are normalized so that perfect prediction has a value of 1, 
%       chance has a value of .5 and having the last prediction be the correct one leads to a value of 0.
%       The .NORMALIZED_RANK_RESULTS has the same .mean_decoding_results and .stdev results 
%       as .ZERO_ONE_LOSS_RESULTS but the different confusion matrix values.  The confusion matrix
%       for the normalized rank results is in the field .confusion_matrix_results.rank_confusion_matrix
%       and contains matrix that is [num_predicted_classes x num_actual_classes x num_training_times x num_test_times]
%       which contains values in the ith row and jth column for the average normalized rank ith predicted class when test points
%       from the jth actual class were presented (i.e., how high up on the predicted labels list is the ith class 
%       when test points from the jth class are shown).  There is also a field .confusion_matrix_results.rank_confusion_matrix_label_mapping 
%       that contains the labels that correspond to the columns of the rank confusion matrix.
%
%    .DECISION_VALUES  This structure contains the decision values.  It has all the same fields as the ZERO_ONE_LOSS_RESULTS
%       except that there is no confusion matrix for this result type.  If save_results.extended_decision_values = 1, then
%       a matrix .classifier_decision_values of dimension [num_resample_runs x num_cv_splits x num_test_points x num_train_times x num_test_times]
%       will be returned that will have the decision values for the correct/actual class for each test point.  
%       If save_results.extended_decision_values = 2, then a matrix .all_classifier_decision_values is returned that has dimensions
%       [num_resample_runs x num_cv_splits x num_classes x num_train_times x num_test_times x num_test_points]
%       that contains all the decision values for every class (not just the correct class).  Thus this result contains 
%       all the information from the whole decoding process (and consequentially it could take up a lot of diskspace/memory).
%       A matrix .all_classifier_decision_labels [num_resample_runs x num_cv_splits x num_test_points x num_train_times x num_test_points] 
%       will also be returned that contains all the labels that were used.  From these two structures it is possible to derive all 
%       other decoding measures that are returned.  
%
%    .ROC_AUC_RESULTS  This structure contains results that measure the area under an receiver operator characteristic (ROC) curves
%       that are created separately for each class from the decision values.  ROC curves graph the proportion of positive test 
%       points correctly classified (true positive rate) as a function of the proportion of negative test points incorrectly
%       classified (false positive rate).  The area under this curve (ROC AUC) gives a measure of decoding accuracy that has a number of
%       useful properties, including the fact that it is invariant to the ratio of positive to negative examples, and that
%       it can be used to determine decoding accuracies when multiple correct classes are present at the same time.  The results
%       in this structure have the following fields:
%   
%       .separate_CV_ROC_results  This structure calculate the ROC AUC separately for each CV split of the data.  The
%           advantage of calculating this separately for each CV split is that it maintains the independence of the
%           decoding results across CV splits.  The disadvantage is that many times there will only be one or a couple of
%           test points for each class which will lead to a highly variable estimate of the ROC AUC.  For this reason it is often
%           better to use the .combined_CV_ROC_results results described below.
%
%       .combined_CV_ROC_results  This structure calculate the ROC AUC by combining all the decision values from all 
%           the test points across the different CV splits.  While the test points should all be independent from one another
%           the classifiers used to evaluate these test points are highly related (since they are most likely trained using very similar
%           data), thus these results could be slightly biased.  However, since a much larger number of points are used to create
%           these ROC curves, the results are more likely to be more sensitive (i.e., these results are slightly more likely to contain
%           type 1 errors but much less likely to create type 2 errors compared to using the .separate_CV_ROC_result structure).
%
%         Both the .separate_CV_ROC_results .combined_CV_ROC_results have the following fields:
%        
%             .decoding_results this contains a  [num_resample_runs x (num_CV_splits) x num_classes x num_training_times x num_test_times] 
%                matrix that contains the ROC AUC results (for the separate_CV_ROC_results there are 5 'dimensions' while for the
%                combined_CV_ROC_results results there are only 4 dimensions since the results combine data from the different CV splits
%                when calculating the ROC AUC results).
%             .mean_decoding_results  This is a [num_training_times x num_test_times] sized matrix that contains the mean ROC AUC
%                values averaged over resample runs the different classes and for the separate_CV_ROC_results, the results are
%                also averaged over CV splits.  
%             .stdev  This field contains the following measure of variability of the ROC AUC results:
%                .over_classes(num_train, num_test, num_resmaples, [num_CV])  Computes the stdev over the results for each class.
%                .over_resamples(num_train, num_test)  Takes the mean over all the classes ([and CV runs]), and computes
%                      the stdev over the resample runs.
%                Additionally, .separate_CV_results has the following measure:
%                  .over_CVs(num_train, num_test, num_resmaples, num_classes)  Takes the mean over classes and computer the variability over the CV runs.
% 
%    .MUTUAL_INFORMATION  This structure contains results that measure the mutual information that is calculated from 
%       the confusion matrix (created from the 0-1 loss results).  For more information on the relationship between
%       mutual information and decoding see Quian Quiroga and Panzeri, Nature Reviews Neuroscience, 2009.  The results
%       in this structure have the following fields:
%
%       .from_combined_confusion_matrix_over_all_resamples  The mutual information in this structure is calculated
%           based on a confusion matrix that combines all the results from all the resample runs (i.e., it is the mutual information 
%           calculated from confusion matrix that is returned in the field ZERO_ONE_LOSS_RESULTS.confusion_matrix_results.confusion_matrix).
%           This mutual information should be the most accurate (and have the least bias during any baseline period), and the accuracy of these
%           results should increase if more resample runs are used.  The disadvantage of this estimate is that no estimate of the variability of this
%           measure is possible since it uses data from all resample runs, thus one cannot plot errorbars with this measure.
%
%       .from_separate_confusion_matrix_for_each_resample  The mutual information in this structure is calculated 
%           based on a confusion matrix that is created separately for each resample run.  This mutual information will most likely be 
%           biased upward due to limited sampling unless a very large number of test points are used (thus the mutual information in the
%           any baseline period is likely to be greater than 0).  The advantage of using this method is that one can now estimate 
%           a measure of this mutual information's variability calculated over resamples, which is in the field .stdev.over_resamples.
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

        datasource;
        classifier;
        feature_preprocessors = [];   % cell array of feature-preprocessor objects
       
        num_resample_runs = 50;
               
        test_only_at_training_times = 0;
        
        
        % result types to save
        save_results = struct('ROC_AUC', 1, 'normalized_rank', 1, 'decision_values', 1, 'extended_decision_values', 0, 'mutual_information', 1);
                              
        % confusion matrix variables
        confusion_matrix_params = struct('create_confusion_matrix', 1, 'create_all_test_points_separate_confusion_matrix', 0, 'save_confusion_matrix_only_train_and_test_at_same_time', 1);
        
        % flags for displaying the progress of the decoding procedure     
        display_progress = struct('resample_run_time', 1, 'zero_one_loss', 1, 'normalized_rank', 0, 'decision_values', 0, 'separate_CV_ROC_results', 0, ...
            'combined_CV_ROC_results', 0, 'training_time_to_display_results', -1, 'convergence_values', 0);    % display how long each resmaple iteration took
                   
        % Setting the different properties of stop_resample_runs_only_when_specfic_results_have_converged causes the resample runs to continue beyond those specified by num_resample runs
        % if the results have not converged to specified values.   
        stop_resample_runs_only_when_specfic_results_have_converged = struct('zero_one_loss_results', [], 'normalized_rank_results', [], 'decision_values', [], 'combined_CV_ROC_results', [], 'separate_CV_ROC_result', [], 'stop_criteria_is_absolute_result_value', 1);
      
                
        
    end
    
    
  
    methods
        
        
        % the constructor  - pass the main objects needed to run the cv decoding
        function cv = standard_resample_CV(the_datasource, the_classifier, the_feature_preprocessors)
                            
            cv.datasource = the_datasource;
            cv.classifier =  the_classifier;
                   
            if nargin > 2    % don't need to specify preprocessing algorithms to run the decoding analysis
                cv.feature_preprocessors = the_feature_preprocessors;   
            end
                      
        end
        
        
 
       
        function DECODING_RESULTS = run_cv_decoding(cv) 

            % sanity checks
            if (cv.confusion_matrix_params.save_confusion_matrix_only_train_and_test_at_same_time == 0) && (cv.test_only_at_training_times == 1)
                error(['To create a confusion matrix that has all training and times one has to run a decoding procedure at all training and test times. '... 
                    'Thus one can not have cv.confusion_matrix_params.save_confusion_matrix_only_train_and_test_at_same_time = 0 and test_only_at_training_times = 1']);
            end
            
            if isempty(cv.test_only_at_training_times)
                cv.test_only_at_training_times = 0;
            end
            

            % might speed things up a bit to assign these to separate variables (matlab is a bit lame)
            datasource = cv.datasource;
            feature_preprocessors = cv.feature_preprocessors;
            classifier = cv.classifier;  
            save_normalized_rank = cv.save_results.normalized_rank;
            save_decision_values = cv.save_results.decision_values;
            save_roc_auc = cv.save_results.ROC_AUC;
            test_only_at_training_times = cv.test_only_at_training_times;
            
            minimum_resample_runs_executed_and_results_have_converged = 0; 
            iResample = 0;
            
            while ~(minimum_resample_runs_executed_and_results_have_converged)   % iResample = 1:cv.num_resample_runs   

                iResample = iResample + 1;
                
                if  cv.display_progress.resample_run_time 
                    tic
                    if iResample == 1, disp('Starting decoding analysis'), end
                end


                % get the data from the datasource
                [all_XTr all_YTr all_XTe all_YTe] = datasource.get_data;  


                % pre-allocating memory for saving stdev.all_single_CV_vals_combined
                if test_only_at_training_times == 1
                    curr_total_CV_zero_one_loss_results = cell(size(all_XTr, 2), 1);
                    curr_total_CV_normalized_rank_results = cell(size(all_XTr, 2), 1);
                    curr_total_CV_decision_values_results = cell(size(all_XTr, 2), 1);                   
                else
                    curr_total_CV_zero_one_loss_results = cell(size(all_XTr, 2), size(all_XTe, 2));
                    curr_total_CV_normalized_rank_results = cell(size(all_XTr, 2), size(all_XTe, 2));
                    curr_total_CV_decision_values_results = cell(size(all_XTr, 2), size(all_XTe, 2));
                end
                
                
                
                
                 % run through all the CV trials first
                for iCV = 1:size(all_XTr{1}, 2)


                    for iTrainingInterval = 1:size(all_XTr, 2)

                            % get the data for the current training interval and CV run
                            XTr = all_XTr{iTrainingInterval}{iCV};
                            
                            if iscell(all_YTr)
                                YTr = all_YTr{iTrainingInterval};
                            else
                                YTr = all_YTr;
                            end
                            
                            
                            % apply preprocessing to the training data 
                            if ~isempty(cv.feature_preprocessors)
                                for iFP = 1:length(cv.feature_preprocessors)
                                    
                                    
                                    [feature_preprocessors{iFP} XTr] = feature_preprocessors{iFP}.set_properties_with_training_data(XTr, YTr);  % save FP parameters and get normalized XTr
                  
                                    
                                    % save preprocessing information (if the user has specified that such information should be saved)

                                      %  [could make the method get_current_info_to_save as not required in the FP interface, 
                                      %    by checking to see if a particular FP has this method via methods(the_feature_preprocessors{iFP}), 
                                      %    but going to keep it a required part of the interface for now]
                                                                                                             
                                    curr_FP_info_to_save = feature_preprocessors{iFP}.get_current_info_to_save;
                                    if ~isempty(curr_FP_info_to_save)
                                        %field_names = fields(curr_FP_info_to_save);  
                                         field_names = fieldnames(curr_FP_info_to_save);  % changed this to make the code compartible with Octave
                                        for iFieldName = 1:length(field_names)
                                            
                                            eval(['curr_FP_data_to_save_one_field = curr_FP_info_to_save.' field_names{iFieldName} ';']); 

                                            if isnumeric(curr_FP_data_to_save_one_field)
                                                eval(['DECODING_RESULTS.FP_INFO{iFP}.'  field_names{iFieldName} '(iResample, iCV, iTrainingInterval, :) = curr_FP_data_to_save_one_field;']);    % will have to modify this slightly if returning a matrix, etc..
                                            else
                                                eval(['DECODING_RESULTS.FP_INFO{iFP}.'  field_names{iFieldName} '{iResample, iCV, iTrainingInterval} = curr_FP_data_to_save_one_field;']); 
                                            end
                                            
                                        end
                                    end

                                     
                                end
                            end   % end preprocessing training code
                                                     
                            
                            
                            % train the classifier 
                            classifier = classifier.train(XTr, YTr);   
                            
                            
                            % if one wants to only test at the same time points that were used for training (to speed things up)
                            if test_only_at_training_times == 1
                                test_interval = iTrainingInterval;
                            else
                                test_interval = 1:size(all_XTe, 2);   % otherwise create full TCT matrix
                            end

                            
                            % run through each test time bin and evaluate how good the decoding accuracy is
                            for iTestInterval = test_interval    
                                
                                if test_only_at_training_times == 1
                                    iTestIntervalSaveInd = 1; 
                                else
                                    iTestIntervalSaveInd = iTestInterval;
                                end
                                
                                
                                % current test data
                                XTe = all_XTe{iTestInterval}{iCV};
                                
                                if iscell(all_YTe)
                                    YTe = all_YTe{iTestInterval};
                                else
                                    YTe = all_YTe;
                                end

                                
                                % apply feature preprocessing to the test data
                                if ~isempty(cv.feature_preprocessors)
                                    for iFP = 1:length(cv.feature_preprocessors)       
                                        XTe = feature_preprocessors{iFP}.preprocess_test_data(XTe);                     
                                    end
                                end


                                % test the classifier
                                [predicted_labels decision_values] = classifier.test(XTe);


                                % store information for creating a confusion matrix
                                if cv.confusion_matrix_params.create_confusion_matrix  
                                    all_predicted_labels_for_confusion_matrix_and_MI(iResample, iCV, :, iTrainingInterval, iTestIntervalSaveInd) = predicted_labels;
                                end

                                
                                % save the zero-one loss results
                                DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.decoding_results(iResample, iCV, iTrainingInterval, iTestIntervalSaveInd, :) = ((length(find(predicted_labels - YTe == 0))/length(predicted_labels)));  
                               
                                
                                % storing stdevs of single 0-1 classifier results
                                DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.all_single_CV_vals(iResample, iCV, iTrainingInterval, iTestIntervalSaveInd) = std(~((predicted_labels - YTe) == 0));   % can also be calculated analytically from decoding_results using p(1-p) binomial variance formula
                                curr_total_CV_zero_one_loss_results{iTrainingInterval,iTestIntervalSaveInd} = [curr_total_CV_zero_one_loss_results{iTrainingInterval,iTestIntervalSaveInd}; ~((predicted_labels - YTe) == 0)];

                                
                                % save the rank_results and/or the decison values...                             
                                if (cv.save_results.decision_values == 1) ||  (cv.save_results.normalized_rank == 1)
                                
                                    % whether one needs to get the rank confusion matrix information (using this variable can slightly speed up the code)
                                    get_confusion_matrix_info = 0;
                                    if (cv.confusion_matrix_params.create_confusion_matrix ~= 0)
                                        if (iTrainingInterval == iTestInterval) 
                                            get_confusion_matrix_info = 1;
                                        else
                                            if ~(cv.confusion_matrix_params.save_confusion_matrix_only_train_and_test_at_same_time == 1), get_confusion_matrix_info = 1; end    
                                        end
                                    end
                                        

                                    % get the normalized rank results, the decision values for the correct class, and rank confusion matrix information                                    
                                    [curr_correct_class_decision_values curr_normalized_rank_results curr_rank_confusion_matrix] = cv.get_rank_and_decision_value_results(YTe, classifier.labels, decision_values, get_confusion_matrix_info);
                                    
                                    
                                   if cv.save_results.normalized_rank == 1
                                    
                                        DECODING_RESULTS.NORMALIZED_RANK_RESULTS.decoding_results(iResample, iCV, iTrainingInterval, iTestIntervalSaveInd) = mean(curr_normalized_rank_results);
                                        DECODING_RESULTS.NORMALIZED_RANK_RESULTS.stdev.all_single_CV_vals(iResample, iCV, iTrainingInterval, iTestIntervalSaveInd) = std(curr_normalized_rank_results);
                                        curr_total_CV_normalized_rank_results{iTrainingInterval,iTestIntervalSaveInd} = [curr_total_CV_normalized_rank_results{iTrainingInterval,iTestIntervalSaveInd}; curr_normalized_rank_results'];   

                                        
                                        % save rank result confusion matrix
                                        if cv.confusion_matrix_params.create_confusion_matrix ~= 0  

                                            YTe_unique_values = unique(YTe);
                                            
                                             % initialize all_rank_results_confusion_matrices
                                             if iCV == 1   
                                                 if  cv.confusion_matrix_params.save_confusion_matrix_only_train_and_test_at_same_time == 1
                                                     curr_resample_rank_confusion_matrix = zeros(length(YTe_unique_values), length(YTe_unique_values), size(all_XTr, 2));
                                                 else
                                                     curr_resample_rank_confusion_matrix = zeros(length(YTe_unique_values), length(YTe_unique_values), size(all_XTr, 2), size(all_XTe, 2));
                                                 end                                                 
                                                 if iResample == 1, total_rank_cm = curr_resample_rank_confusion_matrix; end  % if this is the first run, initialize total_rank_cm matrix 
                                             end

                                            % save rank confusion matrix stuff 
                                            if cv.confusion_matrix_params.save_confusion_matrix_only_train_and_test_at_same_time == 1
                                                if iTrainingInterval == iTestInterval
                                                    curr_resample_rank_confusion_matrix(:, :, iTrainingInterval) = squeeze(curr_resample_rank_confusion_matrix(:, :, iTrainingInterval)) + curr_rank_confusion_matrix;  % using this to get the average predicted ranking for each actual class
                                                end
                                            else
                                                curr_resample_rank_confusion_matrix(:, :, iTrainingInterval, iTestInterval) = squeeze(curr_resample_rank_confusion_matrix(:, :, iTrainingInterval, iTestInterval)) + curr_rank_confusion_matrix;  % using this to get the average predicted ranking for each actual class
                                            end
                                            
                                        end  % end save rank confusion matrix
                                        
                                   end  % end save rank results
                                       
                                    
                                    % save decision value information
                                    if cv.save_results.decision_values == 1
                                        DECODING_RESULTS.DECISION_VALUES.decoding_results(iResample, iCV, iTrainingInterval, iTestIntervalSaveInd) = mean(curr_correct_class_decision_values);  
                                        DECODING_RESULTS.DECISION_VALUES.stdev.all_single_CV_vals(iResample, iCV, iTrainingInterval, iTestIntervalSaveInd) = std(curr_correct_class_decision_values);
                                        curr_total_CV_decision_values_results{iTrainingInterval,iTestIntervalSaveInd} = [curr_total_CV_decision_values_results{iTrainingInterval,iTestIntervalSaveInd}; curr_correct_class_decision_values'];
                                    end
                                    
                                    
                                    if cv.save_results.extended_decision_values == 1
                                        DECODING_RESULTS.DECISION_VALUES.classifier_decision_values(iResample, iCV, :, iTrainingInterval, iTestIntervalSaveInd) = curr_correct_class_decision_values;  
                                    end
                                    
                                    % saving extended decision values for further analysis (this could take up a lot of memory/disk space) 
                                    if cv.save_results.extended_decision_values == 2
                                        DECODING_RESULTS.DECISION_VALUES.all_classifier_decision_values(iResample, iCV, :, :, iTrainingInterval, iTestIntervalSaveInd) = decision_values;  
                                        DECODING_RESULTS.DECISION_VALUES.all_classifier_decision_labels(iResample, iCV, :, iTrainingInterval, iTestIntervalSaveInd) = YTe;  
                                    end

                                    
                                end  % end for saving rank results and/or decision values


                                

                                % save area under ROC curve results
                                if (cv.save_results.ROC_AUC == 1) || (cv.save_results.ROC_AUC == 2) || (cv.save_results.ROC_AUC == 3)

                                    unique_labels = unique(YTe); 
                                    for iClass = 1:size(decision_values, 2)

                                        decision_vals_label_present = decision_values((classifier.labels(iClass) == YTe), iClass);
                                        decision_vals_label_absent = decision_values((classifier.labels(iClass) ~= YTe), iClass);

                                        if (cv.save_results.ROC_AUC == 1) || (cv.save_results.ROC_AUC == 3)  % if saving separate CV ROC results
                                            DECODING_RESULTS.ROC_AUC_RESULTS.separate_CV_ROC_results.decoding_results(iResample, iCV, iClass, iTrainingInterval, iTestIntervalSaveInd) = get_AUC(decision_vals_label_present,  decision_vals_label_absent);
                                        end
                                        
                                        if (cv.save_results.ROC_AUC == 1) || (cv.save_results.ROC_AUC == 2)  % if saving combined CV ROC results 
                                            if iCV == 1
                                                all_CV_decision_vals_label_present{iClass, iTrainingInterval, iTestIntervalSaveInd} = []; 
                                                all_CV_decision_vals_label_absent{iClass, iTrainingInterval, iTestIntervalSaveInd}  = []; 
                                            end
                                            
                                            all_CV_decision_vals_label_present{iClass, iTrainingInterval, iTestIntervalSaveInd}  = [all_CV_decision_vals_label_present{iClass, iTrainingInterval, iTestIntervalSaveInd};  decision_vals_label_present];
                                            all_CV_decision_vals_label_absent{iClass, iTrainingInterval, iTestIntervalSaveInd}  = [all_CV_decision_vals_label_absent{iClass, iTrainingInterval, iTestIntervalSaveInd};  decision_vals_label_absent]; 
                                        end
                                        
                                    end

                                    % save stdev over classes
                                    if (cv.save_results.ROC_AUC == 1) || (cv.save_results.ROC_AUC == 3)
                                        DECODING_RESULTS.ROC_AUC_RESULTS.separate_CV_ROC_results.stdev.over_classes(iResample, iCV, iTrainingInterval, iTestIntervalSaveInd) = std(squeeze(DECODING_RESULTS.ROC_AUC_RESULTS.separate_CV_ROC_results.decoding_results(iResample, iCV, :, iTrainingInterval, iTestIntervalSaveInd)));
                                    end
                                    
                                end  % end save ROC results


                                
                                

                            end  % end for the test time periods

                            
                            
                    end  % end for the train time periods


                    
                end  % end for the CV blocks
                

                
                
                
                % save more measures of decoding result variability
                for iTrainingInterval = 1:size(all_XTr, 2)
                    
                    if cv.test_only_at_training_times == 1  
                        test_interval = iTrainingInterval; 
                    else 
                        test_interval = 1:size(all_XTe, 2); 
                    end
                    
                    for iTestInterval = test_interval 
                                                                                
                        if test_only_at_training_times == 1
                            iTestIntervalSaveInd = 1; 
                        else
                            iTestIntervalSaveInd = iTestInterval;
                        end
                                 
                        DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.all_single_CV_vals_combined(iResample, iTrainingInterval, iTestIntervalSaveInd) = std(curr_total_CV_zero_one_loss_results{iTrainingInterval, iTestIntervalSaveInd});
                        DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_CVs(iResample, iTrainingInterval, iTestIntervalSaveInd) = std(squeeze(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.decoding_results(iResample, :, iTrainingInterval, iTestIntervalSaveInd)));
               
                         if save_normalized_rank == 1
                             DECODING_RESULTS.NORMALIZED_RANK_RESULTS.stdev.all_single_CV_vals_combined(iResample, iTrainingInterval, iTestIntervalSaveInd) = std(curr_total_CV_normalized_rank_results{iTrainingInterval, iTestIntervalSaveInd});
                             DECODING_RESULTS.NORMALIZED_RANK_RESULTS.stdev.over_CVs(iResample, iTrainingInterval, iTestIntervalSaveInd) = std(squeeze(DECODING_RESULTS.NORMALIZED_RANK_RESULTS.decoding_results(iResample, :, iTrainingInterval, iTestIntervalSaveInd)));
                         end
                                                   
                        if save_decision_values == 1 
                            DECODING_RESULTS.DECISION_VALUES.stdev.all_single_CV_vals_combined(iResample, iTrainingInterval, iTestIntervalSaveInd) = std(curr_total_CV_decision_values_results{iTrainingInterval, iTestIntervalSaveInd});
                            DECODING_RESULTS.DECISION_VALUES.stdev.over_CVs(iResample, iTrainingInterval, iTestIntervalSaveInd) = std(squeeze(DECODING_RESULTS.DECISION_VALUES.decoding_results(iResample, :, iTrainingInterval, iTestIntervalSaveInd)));
                        end
                       
                        
                        if (save_roc_auc == 1) || (save_roc_auc == 2) || (save_roc_auc == 3)  
                            if (save_roc_auc == 1) || (save_roc_auc == 2)
                                for iClass = 1:size(decision_values, 2)  
                                      DECODING_RESULTS.ROC_AUC_RESULTS.combined_CV_ROC_results.decoding_results(iResample, iClass, iTrainingInterval, iTestIntervalSaveInd) =  get_AUC(all_CV_decision_vals_label_present{iClass, iTrainingInterval, iTestIntervalSaveInd}, all_CV_decision_vals_label_absent{iClass, iTrainingInterval, iTestIntervalSaveInd});     
                                end  
                                % save combined stdev over classes
                                DECODING_RESULTS.ROC_AUC_RESULTS.combined_CV_ROC_results.stdev.over_classes(iResample, iTrainingInterval, iTestIntervalSaveInd)  = std(squeeze(DECODING_RESULTS.ROC_AUC_RESULTS.combined_CV_ROC_results.decoding_results(iResample, :, iTrainingInterval, iTestIntervalSaveInd)));
                            end
                                
                            % save separate stdevs
                            if (save_roc_auc == 1) || (save_roc_auc == 3)
                                DECODING_RESULTS.ROC_AUC_RESULTS.separate_CV_ROC_results.stdev.over_CVs(iResample, iTrainingInterval, iTestIntervalSaveInd) = std(squeeze(mean(DECODING_RESULTS.ROC_AUC_RESULTS.separate_CV_ROC_results.decoding_results(iResample, :, :, iTrainingInterval, iTestIntervalSaveInd), 2)));
                            end
                        end
                        
                    end 
                    
                end  % end saving more measure of decoding variability 
                
                
                
   
                % save some results for creating the rank confusion matrix
                % normalizing the rank confusion matrix so that each column (actual label) contains the real normalized rank values (i.e., values in the range [0 1])
                if (cv.save_results.normalized_rank == 1) && (cv.confusion_matrix_params.create_confusion_matrix == 1) 

                    if cv.confusion_matrix_params.save_confusion_matrix_only_train_and_test_at_same_time == 1
                         for iCMTime = 1:size(all_XTr, 2)
                             total_num_test_examples_from_each_class = squeeze(sum(curr_resample_rank_confusion_matrix(:, :, iCMTime)))./(length(unique(YTe))./2);   %  (length(unique(YTe))./2) is the sum of all the rank values, i.e., if k = 1/num_clases then this is 0 + 1/k + 2/k + ... + 1
                             curr_resample_rank_confusion_matrix(:, :, iCMTime) = squeeze(curr_resample_rank_confusion_matrix(:, :, iCMTime))./repmat(total_num_test_examples_from_each_class, length(unique(YTe)), 1);
                         end     
                         
                    else
                         for iCMTrain = 1:size(all_XTr, 2)
                             for iCMTest = 1:size(all_XTe, 2)   
                                total_num_test_examples_from_each_class = squeeze(sum(curr_resample_rank_confusion_matrix(:, :, iCMTrain, iCMTest)))./(length(unique(YTe))./2);
                                curr_resample_rank_confusion_matrix(:, :, iCMTrain, iCMTest) = squeeze(curr_resample_rank_confusion_matrix(:, :, iCMTrain, iCMTest))./repmat(total_num_test_examples_from_each_class, length(unique(YTe)), 1);
                             end
                         end     
                    end

                     total_rank_cm = total_rank_cm + curr_resample_rank_confusion_matrix;  

                end

                
                % get convergence values for all the decoding results that are being saved
               if iResample > 1 
                    %result_types = fields(DECODING_RESULTS);
                    result_types = fieldnames(DECODING_RESULTS);  % changed to make the code compartible with Octave
                    for iResultType = 1:length(result_types)

                        if ~strcmp(result_types{iResultType}, 'ROC_AUC_RESULTS')
                            eval(['[max_absolute_diff_ith_deleted percent_change_relative_to_max_ith_deleted] = cv.get_convergence_values(squeeze(mean(DECODING_RESULTS.' result_types{iResultType} '.decoding_results, 2)));'])
                            eval(['convergence_values.' lower(result_types{iResultType}) ' = [max_absolute_diff_ith_deleted percent_change_relative_to_max_ith_deleted];'])
                            if strcmp(result_types{iResultType}, 'ZERO_ONE_LOSS_RESULTS'), convergence_values.zero_one_loss_results(1) = convergence_values.zero_one_loss_results(1) .* 100; end    

                        else
                            if isfield(DECODING_RESULTS.ROC_AUC_RESULTS', 'separate_CV_ROC_results')
                                results_for_each_resample_run_averaged_over_CV_splits = (squeeze(mean(DECODING_RESULTS.ROC_AUC_RESULTS.combined_CV_ROC_results.decoding_results, 2))); % AUROC combined CV values
                                [max_absolute_diff_ith_deleted percent_change_relative_to_max_ith_deleted] = cv.get_convergence_values(results_for_each_resample_run_averaged_over_CV_splits);
                                convergence_values.separate_CV_ROC_result = [max_absolute_diff_ith_deleted percent_change_relative_to_max_ith_deleted];        
                            end
                            if isfield(DECODING_RESULTS.ROC_AUC_RESULTS', 'combined_CV_ROC_results')
                                results_for_each_resample_run_averaged_over_CV_splits = (squeeze(mean(mean(DECODING_RESULTS.ROC_AUC_RESULTS.separate_CV_ROC_results.decoding_results, 3), 2))); % AUROC separate CV values
                                [max_absolute_diff_ith_deleted percent_change_relative_to_max_ith_deleted] = cv.get_convergence_values(results_for_each_resample_run_averaged_over_CV_splits);
                                convergence_values.combined_CV_ROC_results = [max_absolute_diff_ith_deleted percent_change_relative_to_max_ith_deleted];        
                            end
                        end

                    end
                
               else
                    convergence_values = []; % set to empty if on first resample run
               end  % end for getting the convergence values

               
               
                 % display progress of decoding procedure...
                 if cv.display_progress.resample_run_time  == 1
                     
                     reample_run_times(iResample) = toc;   % display how long a resample iteration took
                     
                     curr_date_time = clock;  % current year, month, day, hour...
                     date_conversion_vector = [(24 * 60.^2) 60.^2 60 1];

                     estimated_num_sec_to_complete = (reample_run_times(iResample) .* (cv.num_resample_runs - iResample));  % simple estimate of completion time is: time_of_last_run * num_runs_remaining
                                                                 
                     % estimate the number of days, hours, minutes and second when the code will be done running 
                     estimated_num_sec_to_complete_plus_num_secs_elaped_today = estimated_num_sec_to_complete + curr_date_time(4:end) * date_conversion_vector(2:end)';
                     estimated_completion_time_days_hours_mins_secs = floor(estimated_num_sec_to_complete_plus_num_secs_elaped_today./date_conversion_vector);
                     hour_min_sec_conversion_vector = [24 60 60]; 
                     subtraction_num = 0;
                     for iHourMinSec = 2:4
                         subtraction_num = (subtraction_num + estimated_completion_time_days_hours_mins_secs(iHourMinSec - 1)) .* hour_min_sec_conversion_vector(iHourMinSec - 1);
                         estimated_completion_time_days_hours_mins_secs(iHourMinSec) = estimated_completion_time_days_hours_mins_secs(iHourMinSec) - subtraction_num;
                     end
                     
                     
                     s1 = ['Completed ' num2str(iResample) ' of ' num2str(cv.num_resample_runs) ' resample runs.  The last run took ' num2str(round(reample_run_times(iResample))) ' seconds.']; 

                     if exist('daysadd') == 2   % if the daysadd function from the Financial toolbox exists, print the date/time when the code will be done running 
  
                         % calculate the date in the future when the code will be done running
                         completion_date_string = datestr(daysadd(datenum(curr_date_time(1:3)), estimated_completion_time_days_hours_mins_secs(1)));  % add days remaining to current date to get completion date

                         % print the results to look nice
                         completion_time_string = num2str(mod(estimated_completion_time_days_hours_mins_secs(2), 12));
                         if rem(estimated_completion_time_days_hours_mins_secs(2), 12), suffix = 'PM'; else suffix = 'AM'; end
                         completion_time_string = [completion_time_string ':' num2str(estimated_completion_time_days_hours_mins_secs(3), '%02d') ' ' suffix]; %':' num2str(estimated_completion_time_days_hours_mins_secs(4))];

                         s2 = ['Estimated completition time is for minimum resample runs: ' completion_time_string ' ' completion_date_string ' (the current time is: ' datestr(now, 'HH:MM PM') ')'];                     

                         fprintf('\n%s\n%s\n', s1, s2)
                     
                     else                         
                         fprintf('\n%s\n', s1)                         
                     end
                     
                 end
          
  
                 
    
                 
                 % check if the criteria are met to stop doing resample runs (i.e, at least num_resample runs have been completed and the results have converged)
                 
                 % stop running the resample loop if the minimum number of resample runs has been completed and the results have converged.
                 if iResample < 2   
                     
                 else  % otherwise, see if the results have converged...
                       
                     convergence_display_string = 'convergence values:  ';
                     
                     minimum_resample_runs_executed_and_results_have_converged = 1; 
                     
                     convergence_result_type_names = {'zero_one_loss_results', 'normalized_rank_results', 'decision_values', 'combined_CV_ROC_results', 'separate_CV_ROC_result'};
                     
                     % determine whether to use the absolute convergence threshold or a percentage convergence threshould
                     index_value = 1; thresh_string = ''; 
                     if cv.stop_resample_runs_only_when_specfic_results_have_converged.stop_criteria_is_absolute_result_value == 0, index_value = 2; thresh_string = '%'; end

                     for iConvergenceType = 1:length(convergence_result_type_names)
                        
                         eval(['curr_convergence_values = convergence_values. ' convergence_result_type_names{iConvergenceType} ';'])

                         % check to see if stopping based on convergence values...
                         eval(['curr_convergence_threshold = cv.stop_resample_runs_only_when_specfic_results_have_converged.' convergence_result_type_names{iConvergenceType} ';']);
                         if ~isempty(curr_convergence_threshold)   
                             % if any convergence value is above the threshold specified keep going...   
                             if (curr_convergence_values(index_value) > curr_convergence_threshold)
                                minimum_resample_runs_executed_and_results_have_converged = 0;
                             end                             
                         end
                         
                         % a string to display the convergence results  
                         convergence_display_string = [convergence_display_string convergence_result_type_names{iConvergenceType} ' [value = ' num2str(curr_convergence_values(index_value)) thresh_string  '  threshold = ' num2str(curr_convergence_threshold) thresh_string ']    '];
   
                     end
                     
                     % keep running the decoding loop if the minimum number of resample runs has not yet been executed...                    
                     if (iResample < cv.num_resample_runs)
                        minimum_resample_runs_executed_and_results_have_converged = 0;
                     end
                     
                     
                     
                     
                     % display convergence properies
                     if (cv.display_progress.convergence_values == 1) && iResample > 1                                             
                         fprintf('\n%s\n\n', convergence_display_string)   
                     end

                     
                 end   % end for checking if the results have converged 
                 

                 
                 % display the progress of the results (made this a separate function to make code easier to read  - hopefully this will not slow things down too much passing DECODING_RESULTS as an argument)
                 if (cv.display_progress.zero_one_loss + cv.display_progress.normalized_rank + cv.display_progress.decision_values + cv.display_progress.separate_CV_ROC_results + cv.display_progress.combined_CV_ROC_results) > 0
                     
                     if (cv.display_progress.zero_one_loss == 1) && isfield(DECODING_RESULTS, 'ZERO_ONE_LOSS_RESULTS'),  results_to_display{1} = DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.decoding_results .* 100;  end    % could put an error message if someone tries to display a result type they are not saving, but I will not bother
                     if cv.display_progress.normalized_rank == 1 && isfield(DECODING_RESULTS, 'NORMALIZED_RANK_RESULTS'), results_to_display{2} = DECODING_RESULTS.NORMALIZED_RANK_RESULTS.decoding_results; end
                     if cv.display_progress.decision_values == 1 && isfield(DECODING_RESULTS, 'DECISION_VALUES.decoding_results'),  results_to_display{3} = DECODING_RESULTS.DECISION_VALUES.decoding_results; end
                     if cv.display_progress.separate_CV_ROC_results && isfield(DECODING_RESULTS.ROC_AUC_RESULTS', 'separate_CV_ROC_results'), results_to_display{4} = DECODING_RESULTS.ROC_AUC_RESULTS.separate_CV_ROC_results.decoding_resultss; end
                     if cv.display_progress.combined_CV_ROC_results && isfield(DECODING_RESULTS.ROC_AUC_RESULTS', 'combined_CV_ROC_results'), results_to_display{5} = DECODING_RESULTS.ROC_AUC_RESULTS.combined_CV_ROC_results.decoding_results; end 
                     cv.display_result_progress(results_to_display, cv.display_progress);

                 end
                 
                 
                
             end  % end for the current resample data run


            
            % save additional parameters about the decoding experiment
            % assuming length(YTe) etc. is always the same for all times (which I think has to be true)
            DECODING_RESULTS.CV_PARAMETERS.num_test_points_on_each_CV_run = length(YTe);     
            DECODING_RESULTS.CV_PARAMETERS.num_training_points_on_each_CV_run = length(YTr);
            DECODING_RESULTS.CV_PARAMETERS.dimension_of_data_points = size(XTr, 1); 
            DECODING_RESULTS.CV_PARAMETERS.unique_labels = unique(YTe);
            DECODING_RESULTS.CV_PARAMETERS.num_unique_labels = length(unique(YTe));
            DECODING_RESULTS.CV_PARAMETERS.num_CV_splits = size(all_XTr{1}, 2);
            DECODING_RESULTS.CV_PARAMETERS.num_resample_runs = cv.num_resample_runs;  
            DECODING_RESULTS.CV_PARAMETERS.toolbox_version_number = get_ndt_version;            
            DECODING_RESULTS.CV_PARAMETERS.classifier_name = class(classifier);
            DECODING_RESULTS.CV_PARAMETERS.stop_resample_runs_only_when_specfic_results_have_converged = cv.stop_resample_runs_only_when_specfic_results_have_converged;
            DECODING_RESULTS.CV_PARAMETERS.num_resample_runs_actually_run = iResample;
            DECODING_RESULTS.CV_PARAMETERS.convergence_values = convergence_values;  % might be more appropriate to save this elsewhere, but ok for now...
            for iFP = 1:length(feature_preprocessors), DECODING_RESULTS.CV_PARAMETERS.feature_preprocessor_names{iFP} = class(feature_preprocessors{iFP}); end
            
            % if the datasource has a method get_DS_properties, save the datasource properties returned my this method as well
            % if ismethod(cv.datasource, 'get_DS_properties')
            %    DECODING_RESULTS.DS_PARAMETERS = cv.datasource.get_DS_properties;
            % end
            
            % changed this lines to a try catch statement instead of using the ismethod funciton since the ismethod function does not work properly in Octave
            try
                DECODING_RESULTS.DS_PARAMETERS = cv.datasource.get_DS_properties;
            catch 
            end 
               
            
            
            
            % save additional results (could have this code in the body of this method, but trying to make things more readable)
            DECODING_RESULTS = cv.save_more_decoding_measures(DECODING_RESULTS);
            
            
            % create and confusion matrices and mutual information results                                       
            if (cv.confusion_matrix_params.create_confusion_matrix == 1) || (cv.confusion_matrix_params.create_all_test_points_separate_confusion_matrix == 1) || cv.save_results.mutual_information == 1

                % Using a separate function to create the confusion matrix (this could make the code more memory intensive, but the code is now easier to read).
                % I am assuming YTe is the same for all time periods, which should always be met as far as I can tell
                % confusion matrix => predicted_label x real_label x num_time_periods x num_time_periods
               DECODING_RESULTS = create_confusion_matrices_and_MI(cv, YTe, all_predicted_labels_for_confusion_matrix_and_MI, DECODING_RESULTS);

               if (cv.save_results.normalized_rank == 1) 
                   DECODING_RESULTS.NORMALIZED_RANK_RESULTS.confusion_matrix_results.rank_confusion_matrix = total_rank_cm./cv.num_resample_runs;
                   DECODING_RESULTS.NORMALIZED_RANK_RESULTS.confusion_matrix_results.rank_confusion_matrix_label_mapping = YTe;
               end
               
            end
  
            

        end    % end run_cv_decoding
    
        
 
    end   % end public methods 
    
    
    % should be a private method but for some reason Octave doesn't like when this is private (though Matlab is fine if it's private) so making it a public method :(
    methods (Access = 'public') 
            [maximum_absolute_difference_between_mean_and_ith_deleted_result maximum_percentage_change_between_mean_and_ith_deleted_result] = get_convergence_values(cv, results_for_each_resample_run_averaged_over_CV_splits);        
    end
    
    % few private methods used by the run_cv_decoding method (have these in separate files to make the code slighly easier to read)
    methods (Access = 'private') 
    
        [correct_class_decision_values normalized_rank_results rank_confusion_matrix] = get_rank_and_decision_value_results(cv, YTe, classifier_labels, decision_values, create_rank_confusion_matrix);
        display_result_progress(cv, results_to_display, display_progress);
        DECODING_RESULTS = create_confusion_matrices_and_MI(cv, YTe, all_predicted_labels_for_confusion_matrix_and_MI, DECODING_RESULTS);
        DECODING_RESULTS = save_more_decoding_measures(cv, DECODING_RESULTS);
        
        % for some reason Octave doesn't like this as a private method so I am a public method - but it should never be called outside of this object
        %[maximum_absolute_difference_between_mean_and_ith_deleted_result maximum_percentage_change_between_mean_and_ith_deleted_result] = get_convergence_values(cv, results_for_each_resample_run_averaged_over_CV_splits);
        
    end  % end private methods
        
    
end   % end class

