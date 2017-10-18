classdef plot_standard_results_TCT_object

%  This object allows one to display the results for training and testing at different time points in the experiment 
%   i.e., temporal cross-training (TCT) results, for results that were saved using the standard_resample_CV.run_cv_decoding method.  
%   The results are plotted in an imagesc matrix where the y-axis indicates the time the classifier was trained
%   and the x-axis indicates the time when the classifier was tested.  Additionally, the results can be displayed
%   as a movie, with each frame showing the results for training the classifier at one particular time (indicated 
%   by an line interval at the bottom of the screen), and testing the classifier at all other times (also the results
%   for training and testing at the same time as displayed as a black curve on the figure for reference).
%
%  The constructor for this method takes a string that has the name of the saved results that should be plotted, i.e., 
%   plot_obj = plot_standard_results_TCT_object('result_file_name').  Alternatively, the constructor can take a
%   [num_training_times x num_test_times] matrix with precompiled decoding results and will display a TCT plot with these precompiled results.
%   This can be useful, for example, if one wants to average several results together from different decoding experiments and then display
%   the average results as a TCT plot, etc. 
%
%  The main method of this object is min_and_max_data_range = plot_obj.plot_results which will create a figure that will display the result matrix, 
%   and also optionally a movie showing the results when training the classifier at different points in time. The method also returns
%   the minimum and maxium value of that data which is useful when one wants to plot data on the same range (see the color_result_range property below).  
%   There are several optional parameters that can be set prior to calling the plot_results method which will change how the results
%   are displayed.  These optional properties are:
%   
%    result_type_to_plot (default = 1).  Specifies which result type should be plotted.
%       If this is set to 1, then 0-1 loss results are plotted.
%       If this is set to 2, normalized rank results are plotted. 
%       If this is set to 3, the mean decision values are plotted. 
%       If this is set to 4, ROC_AUC results run separately on each CV split are plotted.    
%       If this is set to 5, ROC_AUC results combined over CV splits are plotted.    
%       If this is set to 6, mutual information created from a confusion matrix that combining data from all resamples is plotted.
%       If this is set to 7, mutual information created from a confusion matrix that is calculated separate and then averaged over resamples is plotted.  
%
%    plot_time_intervals (default []). This property specifies which time points in the experiment the results correspond to 
%      (i.e., this specifies the x-axis values that the results are plotting against).  If this field is empty
%      then this object will attempt to create the x-axis based on the binning_parameter properties that were created by the create_binned_data_from_raster_data       
%      function that is passed through the DS and CV objects (if the binning _parameters structure does not exist, the results will be plotted against sequential
%      integers). If this value is set to a vector, then the times listed on the sides of the TCT matrix will use the values in this time range. 
%      Alternatively, this property can be set to a structure (or a cell array of structs with the same length as the number of results) with the following fields: 
%      plot_time_intervals.bin_width, plot_time_intervals.sampling_interval, and optionally plot_time_intervals.alignment_event_time,
%      and plot_time_intervals.alignment_type,  which will create a time interval that has the corresponding bin widths, step sizes
%      and zero time.  plot_time_intervals.alignment_type specifies whether the bins should be aligned to the center (lot_time_intervals.alignment_type = 2, default value), 
%      the beginning of the bin (plot_time_intervals.alignment_type = 1), or the end of the bin (lot_time_intervals.alignment_type = 3). 
%
%   significant_event_times.  This will cause vertical lines to be drawn at the times specified in this vector which can be used to
%     indicate significant events that occurred during a trial.
%
%   significant_training_event_times. This will cause horizontal lines to be drawn at the times specified in this vector which can be used to
%     indicate significant events that occurred during at particular training times in a trial.
%
%   plot_inds.  Allows plot a smaller range of the results given by the indicies listed in this vector.
%
%   saved_results_structure_name.  By default this object assumes that all results were saved in a structure called 'DECODING_RESULTS'. If
%     the structure of the saved results has a different name, the name can be specified here and the results will be plotted using this name.
%
%   color_result_range.  Specifies the range of the colormap for showing worst to best results.
%
%   figure_position (default = [])  Specifies the position of where the figure show be plotted.  This can be useful because resizing the figure after
%    it is plotted can potentially cause errors in the tick labels.  Also, if a figure is opened to a particular size prior to running this function
%    and this property is left empty, the function will use the size of the current figure.
%
%   TCT_figure_number (default = 1).  The figure number to display the TCT matrix.
%
%   plot_colorbar (default = 1).  If this is set to zero, the colorbar will not be plotted.
%      
%   decoding_result_type_name (default = []).  Sets the name of result type chosen.  If this is left empty, the name will be chosen based on the type
%     of result_type used.  If displaying a movie of the results, this property will also be used for the label y-axis of the movie.
%
%   font_size (default = 16).  The font size for the axis labels.                       
%
%   ylabel_name (default = 'Train time (ms)'.  The label of the y-axis.
%
%   xlabel_name (default = 'Test Time (ms)' ). The label of the x-axis.    
%
%   plot_training_latencies_increasing_up_the_y_axis (default = 1).  If this is set to 1, then the training latencies that are earlier in the trial
%     will be plotted closer to the x-axis (i.e., the results will be equivalent to using 'axis xy' rather than 'axis ij')
%
%
%
%  To display movies of the TCT results, where each frame shows the results for training at t1 and testing at t2, the display_TCT_movie must be set to one 
%   (which is the default behavior of this property).  If plot_obj.display_TCT_movie = 1; then the following additional properties can be set:
%
%  the display_TCT_movie (default = 1).  Plots the temporal-cross-training results as a movie. 
%
%  display_movie_pause_time (default = -1).  How long each frame of the movie show be displayed for.  If this is set to a value less than one, 
%    then the frame will be shown until a key is pressed.
%
%   movie_figure_number (default = 2).  The figure number that the TCT movie will be shown in.
%      
%   movie_save_name (default = []).  If this property is set to a string, the movie will be saved as an avi file using this name.
%
%   movie_time_period_titles  This field allows the movie to display different titles at different time points (which is useful if 
%     one wants to give information about particular events such as 'stimulus on', 'stimulus off', etc.). If this field exist and 
%     if there subfields 'movie_time_period_titles.title_start_times' that lists an array of 
%     times specifying when particular figure titles should be displayed, and and 'movie_time_period_titles.title_names' 
%     that is a cell array that lists the names of the titles at particular times, then the title of the movie figure will 
%     change to display these titles at the given times.  Note that movie_time_period_titles.title_start_times and 
%     movie_time_period_titles.title_names must be the same length.
%
%   chance_level.  Draws a horizontal line at what the chance decoding level is.  If this is unspecified, then the change level will be
%     1/num_classes for 0-1 loss results, no line will be drawn for decision values, and .5 line will be drawn for all other result types.
%
%   line_width (default = 2).  The line width of the plotted movie results.
%
%   sliding_result_color (default = 'b').  The color of the TCT results.  This color is also used
%    for indicating which interval the training time is at the bottom of the plot.  
%
%    errorbar_file_name  If this is set to a string that contains decoding results in 'standard format',
%      errorbars for the movie will be plotted based on the file name listed.  
%
%    errorbar_type_to_plot (default = 1).  If errorbar_file_name is specified, this field specified which type
%      of decoding variance measure (stdev field) should be plotted.  The values for this field are:
%      If this is set to 1, then the standard deviations over resample runs (i.e., stdev.over_resamples) is used.
%      If this is set to 2, then the standard deviations over the mean decoding results from each CV split (i.e., stdev.over_CVs) is used.
%        The mean value over all resample trials is then plotted.
%      If this is set to 3, then the standard deviations over the mean decoding results from each CV split combined together from all resample runs 
%        (i.e., over_CVs_combined_over_resamples) is used.
%      If this is set to 4, then the standard deviation of individual results (e.g., 0 1 values if the 0-1 loss is used) from each CV split (all_single_CV_vals) 
%        are used.  The mean value over all resample runs and CV splits is plotted.
%      If this is set to 5, then the standard deviation of individual results (e.g., 0 1 values if the 0-1 loss is used) from each CV split combined
%        from all CV splits in a single resample run (all_single_CV_vals_combined,) are used.  The mean value over all resample runs is plotted.
%      If this is set to 6, and ROC_AUC results are plotted, then the standard deviaion of the ROC results from over the different classes (over_classes) 
%        are used. The mean results over resample runs (and CV splits for the separate_CV_ROC_results) are plotted.  This value can only be set when 
%        ROC_AUC values are plotted (i.e., result_type_to_plot = 4 or 5)
%      It should also be noted that if separate_CV_ROC_results plotted (i.e., result_type_to_plot = 4), then errorbar_type_to_plot can only be set to values of 1, 2, or 6
%        and if combined_CV_ROC_results plotted (i.e., result_type_to_plot = 5), then errorbar_type_to_plot can only be set to values of 1 or 6.
%
%    errorbar_stdev_multiplication_factor (default = 1).  When plotting the errorbars, the default behavior is to plot then plus and minus 1 stdev for 
%        the stdev type that is specified.  If errorbar_stdev_multiplication_factor is set to a value of k, then the errorbars will be plotted as 
%        mean_results + (k * stdev)  and mean_results - (k * stdev). 
%
%    errorbar_alpha_transparency (default = .2).  Sets the transparency level of the errorbars.
%
%    errorbar_edge_alpha_transparency (default = .2).  Sets the transparency level of the edges of the errorbars.
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


  properties 
      
      result_file_name;
      result_type_to_plot = 1;   % 1: 0-1 loss, 2: normalized rank, 3: mean decision values 4: Area under ROC curve separate, 5: Areau udner ROC curve combined over CV splits
           
      plot_time_intervals;
      significant_event_times = [];
      significant_training_event_times = [];
      plot_inds = [];  % if only want to plot a smaller range of the results 
      color_result_range = [];  
      
      figure_position = [];
      plot_colorbar = 1;
      
      decoding_result_type_name = [];
                       
      saved_results_structure_name = 'DECODING_RESULTS';
      
      TCT_figure_number = 1;

      plot_training_latencies_increasing_up_the_y_axis = 1;
      
      
      % parameters for creating a TCT movie

      display_TCT_movie = 1;
      display_movie_pause_time = -1;

      movie_figure_number = 2;
      
      movie_save_name = [];
      
      movie_time_period_titles = [];
      
      errorbar_file_name = [];
      errorbar_type_to_plot = 1;     % 1: over_resamples, 2: over_CVs, 3: over_CVs_combined_over_resamples, 4: all_single_CV_vals_combined, 5: all_single_CV_vals, 6: over_classes
      errorbar_stdev_multiplication_factor = 1;
      errorbar_alpha_transparency = .2;
      errorbar_edge_alpha_transparency = .2; 
                  
      chance_level = NaN;
      
      sliding_result_color = 'b'; 
      
      line_width = 2;
      font_size = 16;
      ylabel_name = 'Train time (ms)';
      xlabel_name = 'Test time (ms)';     
  
    
  end



    methods 

        % constructor 
        function plot_obj = plot_standard_results_TCT_object(result_file_name)
            
            plot_obj.result_file_name = result_file_name;
            
        end
        

        function min_and_max_data_range = plot_results(plot_obj)
        
            
            if isstr(plot_obj.result_file_name)
            
 
                if plot_obj.result_type_to_plot == 1
                    result_type_name = 'ZERO_ONE_LOSS_RESULTS';
                    if isempty(plot_obj.decoding_result_type_name), plot_obj.decoding_result_type_name = 'Classification Accuracy'; end 
                elseif plot_obj.result_type_to_plot == 2
                    result_type_name = 'NORMALIZED_RANK_RESULTS';
                    if isnan(plot_obj.chance_level), plot_obj.chance_level = .5; end
                    if isempty(plot_obj.decoding_result_type_name), plot_obj.decoding_result_type_name = 'Normalized Rank Accuracy'; end 
                elseif plot_obj.result_type_to_plot == 3
                    result_type_name = 'DECISION_VALUES';
                    if isempty(plot_obj.decoding_result_type_name), plot_obj.decoding_result_type_name = 'Mean Decision Values'; end 
                elseif plot_obj.result_type_to_plot == 4
                    result_type_name = 'ROC_AUC_RESULTS.separate_CV_ROC_results';
                    if isnan(plot_obj.chance_level), plot_obj.chance_level = .5; end
                    if isempty(plot_obj.decoding_result_type_name), plot_obj.decoding_result_type_name = 'Area under ROC curve'; end 
                elseif plot_obj.result_type_to_plot == 5
                    result_type_name = 'ROC_AUC_RESULTS.combined_CV_ROC_results';
                    if isnan(plot_obj.chance_level), plot_obj.chance_level = .5; end
                    if isempty(plot_obj.decoding_result_type_name), plot_obj.decoding_result_type_name = 'Area under ROC curve'; end                 
                elseif plot_obj.result_type_to_plot == 6    
                     result_type_name = 'MUTUAL_INFORMATION.from_combined_confusion_matrix_over_all_resamples';
                     if isnan(plot_obj.chance_level), plot_obj.chance_level = 0; end
                     if isempty(plot_obj.ylabel_name) &&  ~isstr(plot_obj.ylabel_name), plot_obj.ylabel_name = 'Mutual Information (bits)'; end 
                elseif plot_obj.result_type_to_plot == 7    
                     result_type_name = 'MUTUAL_INFORMATION.from_separate_confusion_matrix_for_each_resample';
                     if isnan(plot_obj.chance_level), plot_obj.chance_level = 0; end
                     if isempty(plot_obj.ylabel_name) &&  ~isstr(plot_obj.ylabel_name), plot_obj.ylabel_name = 'Mutual Information (bits)'; end                       
                end

                
                all_results = load(plot_obj.result_file_name);
                all_results_data = eval(['all_results.'  plot_obj.saved_results_structure_name]);

                if plot_obj.result_type_to_plot ~= 6  
                    result_to_plot = eval(['all_results_data.' result_type_name '.mean_decoding_results']);
                else
                    result_to_plot = eval(['all_results_data.' result_type_name '.decoding_results']);
                end


                if size(result_to_plot, 2) == 1
                    error('Can only use this function for results in which the classifier was trained and tested at all times (i.e., .mean_decoding_results must be a matrix');
                end


                % for 0-1 loss results, find chance level if not specified, and also convert results to percent correct from proportion correct
                if (plot_obj.result_type_to_plot == 1)
                    if isnan(plot_obj.chance_level) 
                        plot_obj.chance_level = (1./all_results_data.CV_PARAMETERS.num_unique_labels) .* 100;
                    end   
                    result_to_plot = result_to_plot .* 100;
                end
                
                
                
                % if plot_obj.plot_time_intervals is not set to a vector (or cell array of vectors) of specific times, try to use the binning parameters to
                %    get the time interval to plot the results against
                if isfield(all_results_data, 'DS_PARAMETERS') && isfield(all_results_data.DS_PARAMETERS, 'binned_site_info') && isfield(all_results_data.DS_PARAMETERS.binned_site_info, 'binning_parameters')
                    curr_default_binning_parameters = all_results_data.DS_PARAMETERS.binned_site_info.binning_parameters;
                else
                    curr_default_binning_parameters = [];
                end

                
            elseif ismatrix(plot_obj.result_file_name)
                
                result_to_plot = plot_obj.result_file_name;
                curr_default_binning_parameters = [];
            
            end
        
            

            % if only plotting only a smaller range,
            if ~isempty(plot_obj.plot_inds)
               result_to_plot = result_to_plot(plot_obj.plot_inds, plot_obj.plot_inds);
            end



            the_time_interval  = plot_obj.plot_time_intervals;





            % if plot_obj.plot_time_intervals has been set by the user to a structure with bin_width, sampling_interval, etc.
            %  then create the x-axis time interval from these parameters, or if the_time_interval  is empty and DS_PARAMETERS.binned_site_info.binning_parameters exist,
            %  use these parameters to create a time interval.  Also make sure that plot_time_intervals has not been set to a vector of times.
               %if (isstruct(curr_time_interval) || ~isempty(curr_default_binning_parameters))  &&  ~(ismatrix(curr_time_interval) && ~isempty(curr_time_interval))  
               if (isstruct(the_time_interval) || ~isempty(curr_default_binning_parameters))  ||  ~(ismatrix(the_time_interval) && ~isempty(the_time_interval))  

                   
                time_interval_obj = time_interval_object;
                time_interval_obj.set_parameters_with_binning_stucture(the_time_interval);

                time_interval_obj.set_all_currently_unset_parameters_with_binning_structure(curr_default_binning_parameters);

                the_time_interval = time_interval_obj.get_time_interval;

                if isempty(time_interval_obj.end_time), the_time_interval = the_time_interval(1:length(result_to_plot)); end


            elseif isempty(the_time_interval)  % otherwise if binning_parmaters don't exist (and no interval has been set) just use consecuative numbers on x-axis
                 the_time_interval = 1:length(result_to_plot);   
            end

  

            
            % display the results
            figure(plot_obj.TCT_figure_number)
            hold off
           
            
            if ~isempty(plot_obj.color_result_range)
                imagesc(result_to_plot, plot_obj.color_result_range);
            else
                imagesc(result_to_plot); 
            end
            
            
            
            min_and_max_data_range = [min(min(result_to_plot))  max(max(result_to_plot))];
            
            if plot_obj.plot_colorbar 
                h = colorbar;
                %set(get(h, 'ylabel'), 'String', plot_obj.decoding_result_type_name, 'Rotation', 270, 'VerticalAlignment', 'Bottom', 'FontSize', plot_obj.font_size);
                title(plot_obj.decoding_result_type_name, 'FontSize', plot_obj.font_size);
            end

            
            % put the figure in the correct size before plotting the x and y tick labels since they will get messed up if the figure is resized after this
            if isempty(plot_obj.figure_position),
                plot_obj.figure_position = get(gcf, 'position');
            end           
            set(gcf, 'position', plot_obj.figure_position);

            
            xticks = get(gca, 'XTick');
            xticks(xticks < 1) = [];  xticks(xticks > length(the_time_interval)) = [];  % added to get things to work in Octave
            if min(get(gca, 'XTick')) < 1
                set(gca, 'XTickLabel', [-1 the_time_interval(xticks)]);  % added to get things to work in Octave
            else
                set(gca, 'XTickLabel', the_time_interval(xticks));
            end

            yticks = get(gca, 'YTick');
            yticks(yticks < 1) = [];  yticks(yticks > length(the_time_interval)) = [];  % added to get things to work in Octave
            if min(get(gca, 'YTick')) < 1
                set(gca, 'YTickLabel', [-1 the_time_interval(yticks)]);  % added to get things to work in Octave
            else               
                set(gca, 'YTickLabel', the_time_interval(yticks));            
            end
            
            
            % add lines at significant event times
            significant_event_inds = (((plot_obj.significant_event_times - the_time_interval(1))./(the_time_interval(end) - the_time_interval(1))) .* (size(result_to_plot, 2) - 1)) + 1;
                      
            for iEvent = 1:length(plot_obj.significant_event_times)
               line([significant_event_inds(iEvent), significant_event_inds(iEvent)], get(gca, 'YLim'), 'color', [0 0 0])
            end
            
            
            % add horizontal lines at significant training times
            significant_training_event_inds = (((plot_obj.significant_training_event_times - the_time_interval(1))./(the_time_interval(end) - the_time_interval(1))) .* (size(result_to_plot, 2) - 1)) + 1;
    
            for iEvent = 1:length(plot_obj.significant_training_event_times)
               line(get(gca, 'XLim'), [significant_training_event_inds(iEvent), significant_training_event_inds(iEvent)], 'color', [0 0 0])
            end
            
            ylabel(plot_obj.ylabel_name, 'FontSize', plot_obj.font_size);
            xlabel(plot_obj.xlabel_name, 'FontSize', plot_obj.font_size);
                        
            
            if plot_obj.plot_training_latencies_increasing_up_the_y_axis == 1
                axis xy
            end
            

            %  plot/show a movie of the TCT results
            if plot_obj.display_TCT_movie > 0
            
                
                if ~isempty(plot_obj.movie_save_name)
                   mov = avifile([plot_obj.movie_save_name '.avi'], 'fps', 2); 
                end
                
                
                min_and_max_decoding_vals = [min(min(result_to_plot)), max(max(result_to_plot))];
                
                
                % if plotting errorbars
                if ~isempty(plot_obj.errorbar_file_name)
                
                    % more sanity checks
                    if (plot_obj.errorbar_type_to_plot > 2) && (plot_obj.errorbar_type_to_plot ~=6) && (plot_obj.result_type_to_plot == 4)
                        error('When plotting ROC_AUC_RESULTS.separate_CV_ROC_results one can only plot stdevs over_resamples, over_CVs, and over_classes (i.e., if plot_obj.result_type_to_plot == 4, then plot_obj.errorbar_type_to_plot must be in the range 1, 2, or 6');
                    end

                    if (plot_obj.errorbar_type_to_plot > 1) && (plot_obj.errorbar_type_to_plot ~=6) && (plot_obj.result_type_to_plot == 5)
                        error('When plotting ROC_AUC_RESULTS.combined_CV_ROC_results one can only plot stdevs over_resamples, and over_classes (i.e., if plot_obj.result_type_to_plot == 5, then plot_obj.errorbar_type_to_plot must be 1 or 6');
                    end



                    all_errorbar = load(plot_obj.errorbar_file_name);
                    all_errorbar_data = eval(['all_errorbar.'  plot_obj.saved_results_structure_name]);


                    if  plot_obj.errorbar_type_to_plot == 1
                        stdev_type_name = 'over_resamples';
                        stdev_to_plot = eval(['all_errorbar_data.' result_type_name '.stdev.' stdev_type_name]);
                    elseif  plot_obj.errorbar_type_to_plot == 2
                        stdev_type_name = 'over_CVs';
                        stdev_to_plot = eval(['all_errorbar_data.' result_type_name '.stdev.' stdev_type_name]);
                        stdev_to_plot = squeeze(mean(stdev_to_plot, 1))';
                    elseif  plot_obj.errorbar_type_to_plot == 3
                        stdev_type_name = 'over_CVs_combined_over_resamples';
                        stdev_to_plot = eval(['all_errorbar_data.' result_type_name '.stdev.' stdev_type_name]);
                    elseif  plot_obj.errorbar_type_to_plot == 4
                        stdev_type_name = 'all_single_CV_vals';
                        stdev_to_plot = eval(['all_errorbar_data.' result_type_name '.stdev.' stdev_type_name]);
                        stdev_to_plot = squeeze(mean(mean(stdev_to_plot, 1), 2));
                    elseif  plot_obj.errorbar_type_to_plot == 5
                        stdev_type_name = 'all_single_CV_vals_combined';
                        stdev_to_plot = eval(['all_errorbar_data.' result_type_name '.stdev.' stdev_type_name]);
                        stdev_to_plot = squeeze(mean(stdev_to_plot, 1))';
                    elseif  plot_obj.errorbar_type_to_plot == 6
                        stdev_type_name = 'over_classes';  
                        if plot_obj.result_type_to_plot < 4
                            error('Can only plot stdevs over classes for AUROC results (i.e., if plot_obj.result_type_to_plot is 4 or 5, then plot_obj.errorbar_type_to_plot must be in the range 1-6');
                        end
                        stdev_to_plot = eval(['all_errorbar_data.' result_type_name '.stdev.' stdev_type_name]);  
                        if plot_obj.result_type_to_plot == 5
                            stdev_to_plot = squeeze(mean(stdev_to_plot, 1))';
                        else
                            stdev_to_plot = squeeze(mean(mean(stdev_to_plot, 1), 2));
                        end
                    end


                    if size(stdev_to_plot, 2) == 1
                        error('Can only create errorbars for results in which the classifier was trained and tested at all times (i.e., .mean_decoding_results must be a matrix');
                    end


                    if plot_obj.result_type_to_plot == 1
                        stdev_to_plot = stdev_to_plot .* 100;
                    end


                    % if only plotting only a smaller range,
                    if ~isempty(plot_obj.plot_inds)
                       stdev_to_plot = stdev_to_plot(plot_obj.plot_inds, plot_obj.plot_inds);
                    end


                end
                
                
                
                % for each training time, plot the results at all test times
             
                h = figure(plot_obj.movie_figure_number);
                                
                for iTime = 1:size(result_to_plot, 1)
                

                    hold off            
                    plot(the_time_interval, diag(result_to_plot), 'LineWidth', plot_obj.line_width, 'Color', [0 0 0]); 
                    hold on    
                                        
                    plot(the_time_interval, result_to_plot(iTime, :), 'LineWidth', plot_obj.line_width, 'Color', plot_obj.sliding_result_color);


                    % plot the errorbars if specified 
                    if ~isempty(plot_obj.errorbar_file_name)
                        upper_errorbar = (diag(result_to_plot) + (plot_obj.errorbar_stdev_multiplication_factor .* diag(stdev_to_plot)))';
                        lower_errorbar = (diag(result_to_plot) - (plot_obj.errorbar_stdev_multiplication_factor .* diag(stdev_to_plot)))';

                        fill([the_time_interval fliplr(the_time_interval)],[upper_errorbar fliplr(lower_errorbar)], [0 0 0], 'EdgeColor', [0 0 0], 'EdgeAlpha', plot_obj.errorbar_edge_alpha_transparency, 'FaceAlpha', plot_obj.errorbar_alpha_transparency)

                        upper_errorbar = (result_to_plot(iTime, :) + (plot_obj.errorbar_stdev_multiplication_factor .* stdev_to_plot(iTime, :)));
                        lower_errorbar = (result_to_plot(iTime, :) - (plot_obj.errorbar_stdev_multiplication_factor .* stdev_to_plot(iTime, :)));

                        fill([the_time_interval fliplr(the_time_interval)],[upper_errorbar fliplr(lower_errorbar)], plot_obj.sliding_result_color, 'EdgeColor', plot_obj.sliding_result_color, 'EdgeAlpha', plot_obj.errorbar_edge_alpha_transparency, 'FaceAlpha', plot_obj.errorbar_alpha_transparency)
                    end
                    
                    
                    
                    % plot a spot at the training time...
                    percent_below_min_val_to_plot_spot = .1;  
                    spot_y_position = min_and_max_decoding_vals(1) - (percent_below_min_val_to_plot_spot .* (min_and_max_decoding_vals(2) - min_and_max_decoding_vals(1)));  
                    plot(the_time_interval(iTime), spot_y_position, 'o', 'Color', plot_obj.sliding_result_color, 'MarkerFaceColor', plot_obj.sliding_result_color)
                    
                    if exist('time_interval_obj')
                        bin_width_half_size = round(time_interval_obj.bin_width./2);
                        line([the_time_interval(iTime) - bin_width_half_size, the_time_interval(iTime) + bin_width_half_size], [spot_y_position spot_y_position], 'color', plot_obj.sliding_result_color, 'LineWidth', 2)
                    end
                                        
                    % add axis labels, lines for significant events, etc.
                    xlabel('Time (ms)', 'FontSize', plot_obj.font_size)
                    ylabel(plot_obj.decoding_result_type_name, 'FontSize', plot_obj.font_size)
                    if ~isnan(plot_obj.chance_level)
                        line(get(gca, 'XLim'), [plot_obj.chance_level plot_obj.chance_level], 'color', [0 0 0])
                    end

                    for iEvent = 1:length(plot_obj.significant_event_times)
                       line([plot_obj.significant_event_times(iEvent), plot_obj.significant_event_times(iEvent)], get(gca, 'YLim'), 'color', [0 0 0])
                    end
                    
                    
                    % if different titles are given for differnt time periods, plot the appropriate title
                    if ~isempty (plot_obj.movie_time_period_titles)
                        if ~isfield(plot_obj.movie_time_period_titles, 'title_start_times') ||  ~isfield(plot_obj.movie_time_period_titles, 'title_names') 
                           warning('plot_obj.movie_time_period_titles.title_start_times and plot_obj.movie_time_period_titles.title_names must both be specified in order to plot particular titles at particular times - skipping plotting titles');
                        elseif length(plot_obj.movie_time_period_titles.title_start_times) ~= length(plot_obj.movie_time_period_titles.title_names)
                           warning('plot_obj.movie_time_period_titles.title_start_times and plot_obj.movie_time_period_titles.title_names must both be the same length - skipping plotting titles');
                        else
                            [junk cur_title_ind] =  min(the_time_interval(iTime) - plot_obj.movie_time_period_titles.title_start_times(find((the_time_interval(iTime) - plot_obj.movie_time_period_titles.title_start_times) >= 0)));
                            title(plot_obj.movie_time_period_titles.title_names{cur_title_ind}, 'FontSize', plot_obj.font_size);
                        end
                    end
                    
                    if exist('bin_width_half_size')
                        axis([the_time_interval(1) - bin_width_half_size, the_time_interval(end) + bin_width_half_size, get(gca, 'YLim')])
                    end
                    
                    
                    if ~isempty(plot_obj.movie_save_name)
                            F = getframe(h);  % getframe(gcf);
                            mov = addframe(mov, F);
                    end
                    
                    
                    if plot_obj.display_movie_pause_time < 0
                        pause;  % pause until a key is pressed
                    else
                        pause(plot_obj.display_movie_pause_time);
                    end
                        
                    
                end  % end the for loop over all times


                
                if ~isempty(plot_obj.movie_save_name)
                    mov = close(mov);
                end
                
                
                
            end   % end if for creating the movie
                
                
                
            
            
        end

        
        
        
    end

end   % end classdef








