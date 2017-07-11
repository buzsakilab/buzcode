 


%%   1.  Create strings listing where the toolbox and the tutoral data directories are and add the toolbox to Matlab's path

toolbox_dir = '../../ndt.1.0.0/';
raster_file_directory_name = 'Zhang_Desimone_7objects_raster_data/';

addpath(toolbox_dir)
add_ndt_paths_and_init_rand_generator



%%  2.   Bin the data using one 400 ms bin that starts 100 ms after the stimulus onset


save_prefix_name = 'Binned_Zhang_Desimone_7object_data';
bin_width = 400;
step_size = 400;

start_time = 601;
end_time = 1000;

binned_data_file_name = create_binned_data_from_raster_data(raster_file_directory_name, save_prefix_name, bin_width, step_size, start_time, end_time);





%%  3.  Create a classifier and a feature proprocessor object

the_classifier = max_correlation_coefficient_CL;
the_feature_preprocessors{1} = zscore_normalize_FP;  


%%  4. Let's first train the classifier to discriminate between objects at the upper location, and test the classifier with objects shown at the lower location


%%  4a.  create labels for which exact stimuli (ID plus position) belong in the training set, and which stimuli belong in the test set
id_string_names = {'car', 'couch', 'face', 'kiwi', 'flower', 'guitar', 'hand'};
 
for iID = 1:7   
   the_training_label_names{iID} = {[id_string_names{iID} '_upper']};
   the_test_label_names{iID} = {[id_string_names{iID} '_lower']};
end



%%  4b.  creata a generalization datasource that produces training data at the upper location, and test data at the lower location
num_cv_splits = 18;
 
specific_labels_names_to_use = 'combined_ID_position';  % use the combined ID and position labels
 
ds = generalization_DS(binned_data_file_name, specific_labels_names_to_use, num_cv_splits, the_training_label_names, the_test_label_names);


%%  4c. run a cross-validation decoding analysis that uses the generalization datasource we created to 
%         train a classifier with data from the upper location and test the classifier with data from the lower location

the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);
the_cross_validator.num_resample_runs = 10;

DECODING_RESULTS = the_cross_validator.run_cv_decoding;


% viewing the results suggests that they are above chance  (where chance is .1429)
DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results 





%%  5.  Training and Testing at all locations


mkdir position_invariance_results;  % make a directory to save all the results
num_cv_splits = 18;
 
id_string_names = {'car', 'couch', 'face', 'kiwi', 'flower', 'guitar', 'hand'};
pos_string_names = {'upper', 'middle', 'lower'};
 
for iTrainPosition = 1:3
    
   tic   % print how long it to run the results for training at one position (and testing at all three positions)
    
   for iTestPosition = 1:3
 
      % create the current labels that should be in the training and test sets 
      for iID = 1:7
            the_training_label_names{iID} = {[id_string_names{iID} '_' pos_string_names{iTrainPosition}]};
            the_test_label_names{iID} =  {[id_string_names{iID} '_' pos_string_names{iTestPosition}]};
      end
 
      % create the generalization datasource for training and testing at the current locations
      ds = generalization_DS(binned_data_file_name, specific_labels_names_to_use, num_cv_splits, the_training_label_names, the_test_label_names);       
 
      % create the cross-validator
      the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);
      the_cross_validator.num_resample_runs = 10;
      
      the_cross_validator.display_progress.zero_one_loss = 0;     % let us supress all the output from the cross-validation procedure
      the_cross_validator.display_progress.resample_run_time = 0;
                 
      DECODING_RESULTS = the_cross_validator.run_cv_decoding;    % run the decoding analysis
 
      % save the results
      save_file_name = ['position_invariance_results/Zhang_Desimone_pos_inv_results_train_pos' num2str(iTrainPosition) '_test_pos' num2str(iTestPosition)]; 
      save(save_file_name, 'DECODING_RESULTS')
 
   end
   
   toc
   
end






%% 6.  plot the results


position_names = {'Upper', 'Middle', 'Lower'}


for iTrainPosition = 1:3
    
    
    % load the results from each training and test location
    for iTestPosition = 1:3
        
        load(['position_invariance_results/Zhang_Desimone_pos_inv_results_train_pos' num2str(iTrainPosition) '_test_pos' num2str(iTestPosition)]);
        all_results(iTrainPosition, iTestPosition) = DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results;
    
    end
    
    % create a bar plot for each training lcoation
    subplot(1, 3, iTrainPosition)
    bar(all_results(iTrainPosition, :) .* 100);
    
    title(['Train ' position_names{iTrainPosition}])
    
    ylabel('Classification Accuracy');
    set(gca, 'XTickLabel', position_names)
    xlabel('Test position')
    
    xLims = get(gca, 'XLim');
    line([xLims], [1/7 1/7] .* 100, 'color', [0 0 0]);    % put a line at the chance level of decoding    
    
end


set(gcf, 'position', [247   315   950   300])





