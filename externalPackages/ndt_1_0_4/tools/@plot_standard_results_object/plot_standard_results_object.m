classdef plot_standard_results_object

%  This object allows one to plot the results that are saved when running 
%   the standard_resample_CV.run_cv_decoding method.  The results are plotted
%   as a function of time using results when the classifier was trained and tested
%   at the same time.  
%
%  The constructor for this method takes a cell array that lists the file names
%   of the results that should be plotted, i.e., plot_obj = plot_standard_results_object(result_file_names),
%   where result_file_names = {'file_name1', 'file_name2', ...}.  Alternatively, the constructor can take a [num_results_to_plot x num_time_bins]
%   matrix that has all the decoding results precomputed (if this later option is used, then p_values can not be a cell array 
%   of strings with null distribution directory names, but must be a cell array of vectors with precomputed p-values).  
%
%  The main method of this object is plot_obj.plot_results which will create a figure that will display the results. There 
%   are several optional parameters that can be set prior to calling the plot_results method which will change how the results
%   are displayed.  These optional properties are:
%   
%    result_type_to_plot (default = 1).  Specifies which result type should be plotted.
%       If this is set to 1, then 0-1 loss results are plotted.
%       If this is set to 2, normalized rank results are plotted. 
%       If this is set to 3, the mean decision values are plotted. 
%       If this is set to 4, ROC AUC results run separately on each CV split are plotted.    
%       If this is set to 5, ROC AUC results combined over CV splits are plotted.    
%       If this is set to 6, mutual information created from a confusion matrix that combining data from all resamples is plotted.
%       If this is set to 7, mutual information created from a confusion matrix that is calculated separate and then averaged over resamples is plotted.  
%
%    errorbar_file_names  If this is set to a cell array of string, it will plot
%      errorbars using the data from the file names listed.  The number of file names 
%      in this cell array must be the same as the number of file names in the result_file_names
%      sent to the constructor of this object (i.e., each result to be plotted much have a corresponding
%      errorbar file if errorbars are to be plotted).  Alternatively, errorbar_file_names can be set to a [num_results_to_plot x num_time_bins] sized matrix
%      that contains precomputed errorbars (and the field errorbar_type_to_plot will be ignored). 
%      
%
%    errorbar_type_to_plot (default = 1).  If errorbar_file_names are specified, this field specified which type
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
%        and if combined_CV_ROC_results plotted (i.e., result_type_to_plot = 5), then errorbar_type_to_plot can only be set to values of 1 or 6.  Also if 
%        mutual information created from a confusion matrix that combining data is plotted then one cannot plot errorbars, and if mutual information 
%        created from a confusion matrix that is calculated separate and then averaged over resamples is plotted then errorbar_type_to_plot can only be set to values of 1.
%
%    errorbar_stdev_multiplication_factor (default = 1).  When plotting the errorbars, the default behavior is to plot then plus and minus 1 stdev for 
%        the stdev type that is specified.  If errorbar_stdev_multiplication_factor is set to a value of k, then the errorbars will be plotted as 
%        mean_results + (k * stdev)  and mean_results - (k * stdev). 
%
%    errorbar_transparency_level (default = .2).  Sets the transparency level of the errorbars.
%
%    errorbar_edge_transparency_level (default = .2).  Sets the transparency level of the edges of the errorbars.
%
%    plot_time_intervals (default []). This property specifies which time points in the experiment the results correspond to 
%      (i.e., this specifies the x-axis values that the results are plotting against).  If this field is empty
%      then this object will attempt to create the x-axis based on the binning_parameter properties that were created by the create_binned_data_from_raster_data       
%      function that is passed through the DS and CV objects (if the binning _parameters structure does not exist, the results will be plotted against sequential
%      integers). If this value is set to a vector, then all results specified in result_file_names will use the same time range.  
%      Alternatively, this property can be set to a structure (or a cell array of structs with the same length as the number of results) with the following fields: 
%      plot_time_intervals.bin_width, plot_time_intervals.sampling_interval, and optionally plot_time_intervals.alignment_event_time,
%      and plot_time_intervals.alignment_type,  which will create a time interval that has the corresponding bin widths, step sizes
%      and zero time.  plot_time_intervals.alignment_type specifies whether the bins should be aligned to the center (lot_time_intervals.alignment_type = 2, default value), 
%      the beginning of the bin (plot_time_intervals.alignment_type = 1), or the end of the bin (lot_time_intervals.alignment_type = 3). This property can also cell array 
%      that is the same length as result_file_names, with each entry of the cell array contains a vector of specified times or
%      a structure specifying bin_width, sampling_interval, etc. separately for each result to be plotted.
%
%   plot_inds.  Allows plot a smaller range of the results given by the indicies listed in this vector. This field can also be set to a cell array 
%        that is the same length as result_file_names, with each entry of the cell array containing a vector of indices that should be plotted for each result.
%
%   significant_event_times.  This will cause vertical lines to be drawn at the times specified in this vector which can be used to
%     indicate significant events that occurred during a trial.
%
%   significant_event_regions.  Setting this to a cell array of 2 element vectors (i.e., significant_event_regions{1} = [t1 t2]) will cause 
%     shaded regions to be drawn between t1 and t2.
%
%   significant_event_region_alphas (default = .1).  If significant_event_regions is set to a cell array of 2 dimensional vectors, the shaded regions
%    with have an alpha level of transparency set by this property.  If this property is set to a vector that is the same size as the cell array, then
%    each shaded region can have a different alpha level as determined by the entries in this vector.  
%
%   significant_event_region_colors (default = [0 0 0]).  If significant_event_regions is set to a cell array of 2 dimensional vectors, the shaded regions
%    with have a color specified by the [r g b] vector set here.  If this property is set to a cell array that is the same size as the 
%    significant_event_regions cell array, then each shaded region can have a different color as determined by the 3 dimensional vectors in each cell.
%
%   chance_level.  Draws a horizontal line at what the chance decoding level is.  If this is unspecified, then the change level will be
%     1/num_classes for 0-1 loss results, no line will be drawn for decision values, and .5 line will be drawn for all other result types.
%
%   p_values.  If this is set to a cell array of numbers that is the same length as result_file_names that vectors that of p-values for each time point, then
%     bars will be shown at the bottom of the plot for all times that the p_values are less than p_value_alpha_level.  If this is set to a cell array of strings
%     that contain the names of directories that have files that comprise a null distribution, the p-values will be calculated and a bars will be shown at the 
%     bottom of the plot for all times that the p_values are less than p_value_alpha_level. The width of the interval indicating whether a data point is above chance
%     is the time when a data point is above chance +/- 1/2 of the bin width.  For more information about how to run decoding analyses to estimate 
%     p-values, see the file tools/create_pvalues_from_null_distribution_files.m).
%
%   p_value_alpha_level.  If the p_values field is set, then  bars will be shown at the bottom of the plot for all times that the p_values are less than this value.
%      
%   add_pvalue_latency_to_legends_alignment (default = 4).  If this property is set to value greater than 0, and if the p_values field has been set, then the
%     latency when the decoding results are above chance will be added to the legend.  If this value is set to 1 it will show the beginning of the time interval when 
%     the results are above chance, if it is set to 2 it will show the middle of the time interval when the results are above chance, if this is set to 3 it will
%     show the end of the time interval when the results are above chance and if this is set to 4 (default) it will show the beginning and the end of the 
%     interval when the results are above chance.  If this property is set to 5, then the time interval that was added to make the results above chance will be displayed
%     (i.e., the time of the first interval that is above chance will be displayed minus any time that is in the previous interval (which was the
%     last time interval that was not above chance)).  
%
%   collapse_all_times_when_estimating_pvals (default = 0). If this is set to one and p_values is set to a cell array that contains the names of directories 
%           with shuffled results, then the p-value will be calculated froma  null distribution in which data from  all time bins are combined together.  
%           The advantage of using this is that if the null distributions at each point in time are the same, then one can get a more precise estimate of the p-values 
%           for the same computational cost. [Added with NDT version 1.4] 
%
%   saved_results_structure_name.  By default this object assumes that all results were saved in a structure called 'DECODING_RESULTS'.  If
%     the structure of the saved results has a different name, the name can be specified here and the results will be plotted using this name.
%
%   the_axis (default = []).  This has the same effect as using the axis([x1 x2 y1 y2]) function but setting this first will allow significant time lines and
%     shaded regions to be drawn to span the whole display, as opposed to if this is set later using the axis function.  
%
%   the_colors.  A cell array listing the colors that the different results should be plotted as.  There are 20 colors defined as the default colors,
%     if more than 20 results are being compared on the same figure, this property must be extended to include more colors.
%
%   line_width (default = 2).  The line width of the plotted results.
%
%   line_style (default = '-').  The line style used (i.e., can use dashed lines, circles, etc.).  
%
%   font_size (default = 16).  The font size for the axis labels.
%
%   ylabel_name.  The default behavior is to display the appropriate name for the given result_type that is being plotted, but this can be set
%   to display alternative ylabels.
%
%   xlabel_name (default = 'Time (ms)' ).  The label of the x-axis.      
%  
%   legend_names.  Allows one to set a cell array of strings corresponding to the legends for the different results
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
      
      result_file_names;
      result_type_to_plot = 1;   % 1: 0-1 loss, 2: normalized rank, 3: mean decision values 4: Area under ROC curve separate, 5: Areau udner ROC curve combined over CV splits
      
      errorbar_file_names = [];
      errorbar_type_to_plot = 1; % 1: over_resamples, 2: over_CVs, 3: over_CVs_combined_over_resamples, 4: all_single_CV_vals_combined, 5: all_single_CV_vals, 6: over_classes
      errorbar_stdev_multiplication_factor = 1;
      errorbar_transparency_level = .2;
      errorbar_edge_transparency_level = .2; 
      
      plot_time_intervals;
      plot_inds = [];  % if only want to plot a smaller range of the results 
      significant_event_times = [];
      
      significant_event_regions = [];
      significant_event_region_alphas = .1;
      significant_event_region_colors = [0 0 0];
      
      p_values = [];
      p_value_alpha_level = eps;
      add_pvalue_latency_to_legends_alignment = 4;   % by default add full start and end interval latency to the legend if p-values are specified and fields .bin_width and .sampling_interval are specified for plot_time_intervals
      collapse_all_times_when_estimating_pvals = 0;  
      
      chance_level = NaN;
                 
      saved_results_structure_name = 'DECODING_RESULTS';
      
      the_axis = [];
      
      %the_colors = {[0 0 1], [1 0 0], [0 1 0], [0 1 1], [1 0 1], [1 1 0], [0 0 0]};
      the_colors = {[0.00  0.00  1.00], [1.00  0.00  0.00], [0.00  0.50  0.00], [0.00  0.75  0.75], [0.75  0.00  0.75], [0.75  0.75  0.00], [0.25  0.25  0.25], ...
            [0.75  0.25  0.25], [0.95  0.95  0.00], [0.25  0.25  0.75], [0.75  0.75  0.75], [0.00  1.00  0.00], [0.76  0.57  0.17], [0.54  0.63  0.22], ...
            [0.34  0.57  0.92], [1.00  0.10  0.60], [0.88  0.75  0.73],[0.10  0.49  0.47], [0.66  0.34  0.65], [0.99  0.41  0.23]};
      
      line_width = 2;
      line_style = '-';
      font_size = 16;
      ylabel_name = [];
      xlabel_name = 'Time (ms)';      
  
      legend_names = {};
      
      
      
  end



    methods 

        % constructor 
        function plot_obj = plot_standard_results_object(result_file_names)
            
            plot_obj.result_file_names = result_file_names;
            
        end
        

        
        function plot_results(plot_obj)
            
            % sanity check that plot_obj.result_file_names is a cell array of strings
            if ~(iscell(plot_obj.result_file_names) || ismatrix(plot_obj.result_file_names))
                error(['plot_obj.result_file_names must be set to a cell array - where each element of the cell array is a string that is the name of a saved DECODING_RESULTS file' ...
                    ' - or plot_obj.result_file_names must be a [num_results x num_time_points] matrix that has precompiled decoding results'])
            end

            
            % get the count of how many results there are to be plotted
            if iscell(plot_obj.result_file_names)
                number_of_results_to_plot = length(plot_obj.result_file_names);
            else
                number_of_results_to_plot = size(plot_obj.result_file_names, 1);
            end

            
            % sanity check that there are enough colors
            if number_of_results_to_plot > length(plot_obj.the_colors)
                error('More results being plotted that colors specified - plot_obj.the_colors need to be set to include more colors');
            end
                            
            
            
            if plot_obj.result_type_to_plot == 1
                result_type_name = 'ZERO_ONE_LOSS_RESULTS';
                if isempty(plot_obj.ylabel_name) &&  ~isstr(plot_obj.ylabel_name), plot_obj.ylabel_name = 'Classification Accuracy'; end 
            elseif plot_obj.result_type_to_plot == 2
                result_type_name = 'NORMALIZED_RANK_RESULTS';
                if isnan(plot_obj.chance_level), plot_obj.chance_level = .5; end
                if isempty(plot_obj.ylabel_name) &&  ~isstr(plot_obj.ylabel_name), plot_obj.ylabel_name = 'Normalized Rank Accuracy'; end 
            elseif plot_obj.result_type_to_plot == 3
                result_type_name = 'DECISION_VALUES';
                if isempty(plot_obj.ylabel_name) &&  ~isstr(plot_obj.ylabel_name), plot_obj.ylabel_name = 'Mean Decision Values'; end 
            elseif plot_obj.result_type_to_plot == 4
                result_type_name = 'ROC_AUC_RESULTS.separate_CV_ROC_results';
                if isnan(plot_obj.chance_level), plot_obj.chance_level = .5; end
                if isempty(plot_obj.ylabel_name) &&  ~isstr(plot_obj.ylabel_name), plot_obj.ylabel_name = 'Area under ROC curve'; end 
            elseif plot_obj.result_type_to_plot == 5
                result_type_name = 'ROC_AUC_RESULTS.combined_CV_ROC_results';
                if isnan(plot_obj.chance_level), plot_obj.chance_level = .5; end
                if isempty(plot_obj.ylabel_name) &&  ~isstr(plot_obj.ylabel_name), plot_obj.ylabel_name = 'Area under ROC curve'; end 
            elseif plot_obj.result_type_to_plot == 6    
                 result_type_name = 'MUTUAL_INFORMATION.from_combined_confusion_matrix_over_all_resamples';
                 if isnan(plot_obj.chance_level), plot_obj.chance_level = 0; end
                 if isempty(plot_obj.ylabel_name) &&  ~isstr(plot_obj.ylabel_name), plot_obj.ylabel_name = 'Mutual Information (bits)'; end 
            elseif plot_obj.result_type_to_plot == 7    
                 result_type_name = 'MUTUAL_INFORMATION.from_separate_confusion_matrix_for_each_resample';
                 if isnan(plot_obj.chance_level), plot_obj.chance_level = 0; end
                 if isempty(plot_obj.ylabel_name) &&  ~isstr(plot_obj.ylabel_name), plot_obj.ylabel_name = 'Mutual Information (bits)'; end                      
            end
            
            
            
            for iResult = 1:number_of_results_to_plot
                

                            
                % if plot_obj.result_file_names is a cell array listing the names of results in standard results format, then load these results
                if iscell(plot_obj.result_file_names)
            

                    all_curr_results = load(plot_obj.result_file_names{iResult});
                    all_curr_results_data = eval(['all_curr_results.'  plot_obj.saved_results_structure_name]);

                    if plot_obj.result_type_to_plot ~= 6  
                        curr_result_to_plot = eval(['all_curr_results_data.' result_type_name '.mean_decoding_results']);
                    else
                        curr_result_to_plot = eval(['all_curr_results_data.' result_type_name '.decoding_results']);  % MI using confusion matrix over all bootstraps does not involve any averaging
                    end


                    if size(curr_result_to_plot, 2) ~= 1
                        curr_result_to_plot = diag(curr_result_to_plot);
                    end


                    % for 0-1 loss results, find chance level if not specified, and also convert results to percent correct from proportion correct
                    if (plot_obj.result_type_to_plot == 1)
                        if isnan(plot_obj.chance_level) 
                            plot_obj.chance_level = (1./all_curr_results_data.CV_PARAMETERS.num_unique_labels) .* 100;
                        end   
                        curr_result_to_plot = curr_result_to_plot .* 100;
                    end


                    % if plot_obj.plot_time_intervals is not set to a vector (or cell array of vectors) of specific times, try to use the binning parameters to 
                    %    get the time interval to plot the results against
                    if isfield(all_curr_results_data, 'DS_PARAMETERS') && isfield(all_curr_results_data.DS_PARAMETERS, 'binned_site_info') && ~isempty(all_curr_results_data.DS_PARAMETERS.binned_site_info)
                        curr_default_binning_parameters = all_curr_results_data.DS_PARAMETERS.binned_site_info.binning_parameters;
                        all_default_binning_parameters{iResult} = all_curr_results_data.DS_PARAMETERS.binned_site_info.binning_parameters;
                    else
                        curr_default_binning_parameters = [];
                    end

                
                         
                else
                    
                     % if plot_obj.result_file_names is a matrix of precomputed results, just return the current row of results
                     curr_result_to_plot = plot_obj.result_file_names(iResult, :)';
                     curr_default_binning_parameters = [];
                     if isempty(plot_obj.ylabel_name) &&  ~isstr(plot_obj.ylabel_name), plot_obj.ylabel_name = 'Custom decoding type'; end 

                end
                    
                
                
                
                % get the time (ususually in milliseconds) where the data should be plotted    

                % if a cell array of vectors/structures has been set, get the current one, otherwise use the same one for all times
                if iscell(plot_obj.plot_time_intervals)
                    curr_time_interval = plot_obj.plot_time_intervals{iResult};
                else
                    curr_time_interval = plot_obj.plot_time_intervals;
                end
  

               % if plot_obj.plot_time_intervals has been set by the user to a structure with bin_width, sampling_interval, etc.
               %  then create the x-axis time interval from these parameters, or if curr_time_interval is empty and DS_PARAMETERS.binned_site_info.binning_parameters exist,
               %  use these parameters to create a time interval.  Also make sure that plot_time_intervals has not been set to a vector of times.
               %if (isstruct(curr_time_interval) || ~isempty(curr_default_binning_parameters))  ||  ~(ismatrix(curr_time_interval) && ~isempty(curr_time_interval))  
               %if ~iscell(plot_obj.plot_time_intervals) && ((isstruct(curr_time_interval) || ~isempty(curr_default_binning_parameters))  ||  ~(ismatrix(curr_time_interval) && ~isempty(curr_time_interval)))  
               if ~iscell(plot_obj.plot_time_intervals) && ((isstruct(curr_time_interval) || ~isempty(curr_default_binning_parameters))  ||  ~(ismatrix(curr_time_interval) && isempty(curr_time_interval)))    % hopefully this is right now :)

                   time_interval_obj = time_interval_object;
                    time_interval_obj.set_parameters_with_binning_stucture(curr_time_interval);
                    
                    time_interval_obj.set_all_currently_unset_parameters_with_binning_structure(curr_default_binning_parameters);
                    
                    curr_time_interval = time_interval_obj.get_time_interval;
                    
                    % if end time is not given, reduce the interval to the proper length  (might not actually be necessary to do this, but doing it anyway)
                    if isempty(time_interval_obj.end_time), curr_time_interval = curr_time_interval(1:length(curr_result_to_plot));  end

                    all_time_interval_objects{iResult} = time_interval_obj;
                    

                elseif isempty(curr_time_interval)  % otherwise if binning_parmaters don't exist (and no interval has been set) just use consecuative numbers on x-axis
                     curr_time_interval = 1:length(curr_result_to_plot);   
                end
                    

                % if only plotting only a smaller range of results, get the inds to use
                if ~isempty(plot_obj.plot_inds)
                   if iscell(plot_obj.plot_inds)
                        curr_plot_inds = plot_obj.plot_inds{iResult};
                   else
                        curr_plot_inds = plot_obj.plot_inds;
                   end
                else
                        curr_plot_inds = 1:length(curr_result_to_plot);
                end
               
                
                % plot the results
                if iResult == 1, hold off; end             
                plot(curr_time_interval(curr_plot_inds), curr_result_to_plot(curr_plot_inds), plot_obj.line_style, 'LineWidth', plot_obj.line_width, 'Color', plot_obj.the_colors{iResult}); 
                hold on    
                
                all_results_to_plot{iResult} = curr_result_to_plot;   % keep the full length results, only truncate with inds_to_use when the results are being plotted
                all_time_intervals{iResult}  = curr_time_interval;
                all_plot_inds{iResult} = curr_plot_inds;
                
            end
            
            

            
            % plot errorbar if given ...   
            if ~isempty(plot_obj.errorbar_file_names)
                
                if iscell(plot_obj.errorbar_file_names)
                    number_of_errorbars = length(plot_obj.errorbar_file_names);
                else
                    number_of_errorbars = size(plot_obj.errorbar_file_names, 1);
                end
                               
                % more sanity checks
                if number_of_errorbars ~= number_of_results_to_plot     %size(plot_obj.errorbar_file_names) ~= number_of_results_to_plot
                    error('There must be as many errorbar_file_names as there are result_file_names')
                end
                                
  
                for iErrorBar = 1:number_of_errorbars


                    if iscell(plot_obj.errorbar_file_names)

                        if (plot_obj.errorbar_type_to_plot > 2) && (plot_obj.errorbar_type_to_plot ~=6) && (plot_obj.result_type_to_plot == 4)
                            error('When plotting ROC_AUC_RESULTS.separate_CV_ROC_results one can only plot stdevs over_resamples, over_CVs, and over_classes (i.e., if plot_obj.result_type_to_plot == 4, then plot_obj.errorbar_type_to_plot must be in the range 1, 2, or 6');
                        end

                        if (plot_obj.errorbar_type_to_plot > 1) && (plot_obj.errorbar_type_to_plot ~=6) && (plot_obj.result_type_to_plot == 5)
                            error('When plotting ROC_AUC_RESULTS.combined_CV_ROC_results one can only plot stdevs over_resamples, and over_classes (i.e., if plot_obj.result_type_to_plot == 5, then plot_obj.errorbar_type_to_plot must be 1 or 6');
                        end

                        if (plot_obj.result_type_to_plot == 6)
                            error('No errorbars can be estimated when plotting MUTUAL INFORMAION from a confusion matrix that is created by combining data from all resample runs (i.e., if plot_obj.result_type_to_plot == 6 then no errorbars can be plotted')
                        end

                        if (plot_obj.result_type_to_plot == 7) && (plot_obj.errorbar_type_to_plot > 1)
                            error('When plotting MUTUAL_INFORMATION calculated from confusion matrices that were computed separately for each resample run, one can only plot stdevs over_resamples,  (i.e., if plot_obj.result_type_to_plot == 7, then plot_obj.errorbar_type_to_plot must be 1');
                        end



                        all_curr_errorbar = load(plot_obj.errorbar_file_names{iErrorBar});
                        all_curr_errorbar_data = eval(['all_curr_errorbar.'  plot_obj.saved_results_structure_name]);


                        if  plot_obj.errorbar_type_to_plot == 1
                            stdev_type_name = 'over_resamples';
                            curr_stdev_to_plot = eval(['all_curr_errorbar_data.' result_type_name '.stdev.' stdev_type_name]);
                        elseif  plot_obj.errorbar_type_to_plot == 2
                            stdev_type_name = 'over_CVs';
                            curr_stdev_to_plot = eval(['all_curr_errorbar_data.' result_type_name '.stdev.' stdev_type_name]);
                            curr_stdev_to_plot = squeeze(mean(curr_stdev_to_plot, 1))';
                        elseif  plot_obj.errorbar_type_to_plot == 3
                            stdev_type_name = 'over_CVs_combined_over_resamples';
                            curr_stdev_to_plot = eval(['all_curr_errorbar_data.' result_type_name '.stdev.' stdev_type_name]);
                        elseif  plot_obj.errorbar_type_to_plot == 4
                            stdev_type_name = 'all_single_CV_vals';
                            curr_stdev_to_plot = eval(['all_curr_errorbar_data.' result_type_name '.stdev.' stdev_type_name]);
                            curr_stdev_to_plot = squeeze(mean(mean(curr_stdev_to_plot, 1), 2));
                        elseif  plot_obj.errorbar_type_to_plot == 5
                            stdev_type_name = 'all_single_CV_vals_combined';
                            curr_stdev_to_plot = eval(['all_curr_errorbar_data.' result_type_name '.stdev.' stdev_type_name]);
                            curr_stdev_to_plot = squeeze(mean(curr_stdev_to_plot, 1))';
                        elseif  plot_obj.errorbar_type_to_plot == 6
                            stdev_type_name = 'over_classes';  
                            if plot_obj.result_type_to_plot < 4
                                error('Can only plot stdevs over classes for AUROC results (i.e., if plot_obj.result_type_to_plot is 4 or 5, then plot_obj.errorbar_type_to_plot must be in the range 1-6');
                            end
                            curr_stdev_to_plot = eval(['all_curr_errorbar_data.' result_type_name '.stdev.' stdev_type_name]);  
                            if plot_obj.result_type_to_plot == 5
                                curr_stdev_to_plot = squeeze(mean(curr_stdev_to_plot, 1))';
                            else
                                curr_stdev_to_plot = squeeze(mean(mean(curr_stdev_to_plot, 1), 2));
                            end
                        end


                        if size(curr_stdev_to_plot, 2) ~= 1
                           curr_stdev_to_plot = diag(curr_stdev_to_plot); 
                        end


                        if plot_obj.result_type_to_plot == 1
                            curr_stdev_to_plot = curr_stdev_to_plot .* 100;
                        end

                        
                    else   % if using precomputed errorbars...                        
                        curr_stdev_to_plot = plot_obj.errorbar_file_names(iErrorBar, :)';
                    end
                                              

                    % actually plot the errorbars
                    upper_errorbar = (all_results_to_plot{iErrorBar} + (plot_obj.errorbar_stdev_multiplication_factor .* curr_stdev_to_plot))';
                    lower_errorbar = (all_results_to_plot{iErrorBar} - (plot_obj.errorbar_stdev_multiplication_factor .* curr_stdev_to_plot))';
                    
                    fill([all_time_intervals{iErrorBar}(all_plot_inds{iErrorBar}) fliplr(all_time_intervals{iErrorBar}(all_plot_inds{iErrorBar}))],[upper_errorbar(all_plot_inds{iErrorBar}) fliplr(lower_errorbar(all_plot_inds{iErrorBar}))], plot_obj.the_colors{iErrorBar}, 'EdgeColor', plot_obj.the_colors{iErrorBar}, 'EdgeAlpha', plot_obj.errorbar_edge_transparency_level, 'FaceAlpha', plot_obj.errorbar_transparency_level)

                end

            end   % end for plotted the errorbars
            
            
            
                
            
            % plot significant times if p-values are given
            if ~isempty(plot_obj.p_values)                
                
                ylims = get(gca, 'YLim');
                y_interval_length = ylims(2) - ylims(1);
                
                % sanity check to make sure enough pvalues have been given
                if number_of_results_to_plot ~= length(plot_obj.p_values)
                    error('plot_obj.p_values must be a cell array that is the same number of results as plot_obj.result_file_names')
                end
                
                
                for iResult = 1:number_of_results_to_plot

                    y_offset = ylims(1) - (((iResult -1) .* .025) .* (y_interval_length));
      
                    % if plot_obj.p_values is a directory name with the null distribution files, calculate p-values from the null distribution
                    if isstr(plot_obj.p_values{iResult})
                        
                        
                      % give an error if one tries to compute p-values from a null distribution file names when customized precomputed decoding results are given
                      %  (could possibly change this so it is still possible to compute pvalues from a directory of null distribution files by changing the p-value object
                      %   (and also trusting that the user knows what he/she is doing) but not going to do this right now).
                       if ~iscell(plot_obj.result_file_names)   
                            error(['If plot_obj.result_file_names is a [num_results x num_time_bins] matrix with precomputed decoding accuracies ' ...
                              'then one can not compute vaules by having plot_obj.p_values{iResult} contain the name of a directory']); 
                       end
                        
                        
                       pval_obj = pvalue_object(plot_obj.result_file_names{iResult}, plot_obj.p_values{iResult});
                       pval_obj.the_result_type = plot_obj.result_type_to_plot;
                       pval_obj.collapse_all_times_when_estimating_pvals = plot_obj.collapse_all_times_when_estimating_pvals;
                       plot_obj.p_values{iResult} = pval_obj.create_pvalues_from_nulldist_files;
                       
                    else   % otherwise just use the pvalues given
                        
                       pval_obj = pvalue_object;
                       pval_obj.p_values = plot_obj.p_values{iResult};
                        
                    end
                    
                    
                    % set the val_obj.p_values to the plot_obj pvalues in case pvalues were passed to the plot_object instead of calculated by the pvalue_object
                    pval_obj.p_values = plot_obj.p_values{iResult};
                    
                    
                   % add decoded information latency to legends
                   if plot_obj.add_pvalue_latency_to_legends_alignment > 0

                       pval_obj.latency_alpha_significance_level = plot_obj.p_value_alpha_level;

                       % could set the pval_obj.latency_num_consecutive_sig_bins here


                       if ~exist('all_time_interval_objects') ||  isempty(all_time_interval_objects{iResult})
                               pval_obj.latency_time_interval = all_time_intervals{iResult};
                       else
                               pval_obj.latency_time_interval = all_time_interval_objects{iResult};
                               pval_obj.latency_time_interval_alignment_type = plot_obj.add_pvalue_latency_to_legends_alignment;
                       end

                       latency_time = pval_obj.get_latency_of_decoded_information;

                       if isscalar(latency_time)
                           curr_pval_latency_string = num2str(latency_time);
                       else
                           curr_pval_latency_string = [num2str(latency_time(1)) '-' num2str(latency_time(2))];
                       end


                       % add the latency string to the legend...
                       if isempty(plot_obj.legend_names) || (numel(plot_obj.legend_names) < iResult)
                           plot_obj.legend_names{iResult} = [ '(' curr_pval_latency_string ')'];
                       else
                           plot_obj.legend_names{iResult} = [plot_obj.legend_names{iResult} ' (' curr_pval_latency_string ')'];
                       end
                           

                   end   % end for adding decoding information latency to the legends

                           

                    % The length of p-value significance markers are based on the distance between adjacent sampled points (i.e., the sampling_interval).  The 
                    %   length of each marker is from the center of the bin to 1/2 the sampling_interval to the left and 1/2 the sampling_interval to the right. 
                    edges = diff(all_time_intervals{iResult}(all_plot_inds{iResult}))./2;
                    start_intervals = all_time_intervals{iResult}(all_plot_inds{iResult}) -  [edges(1) edges];   % for first interval, left marker length be equal to the right marker length
                    end_intervals = all_time_intervals{iResult}(all_plot_inds{iResult}) + [edges edges(end)];   % for the last interval, have the right marker length be equal to the left marker length
                    inds_to_use = find(plot_obj.p_values{iResult}(all_plot_inds{iResult}) < plot_obj.p_value_alpha_level); 
                                      
                    for iSigInterval = 1:length(inds_to_use)
                        line([start_intervals(inds_to_use(iSigInterval)) end_intervals(inds_to_use(iSigInterval))], [y_offset y_offset], 'color', plot_obj.the_colors{iResult}, 'LineWidth', 5)      
                    end      
                    
                    
                end
                
                axis([get(gca, 'XLim') (ylims(1) - (((iResult + 1.5) .* .025) .* (y_interval_length))) ylims(2)])
                
            end
              
            
            % add axis labels, lines for significant events, etc.
            xlabel(plot_obj.xlabel_name, 'FontSize', plot_obj.font_size)
            ylabel(plot_obj.ylabel_name, 'FontSize', plot_obj.font_size)
            if ~isnan(plot_obj.chance_level)
                line(get(gca, 'XLim'), [plot_obj.chance_level plot_obj.chance_level], 'color', [0 0 0])
            end
            
            
            % if axes length have been specified, make the figure have the appropriate axes before plotting significant even times/regions
            if ~isempty(plot_obj.the_axis)
                axis(plot_obj.the_axis);  
            end
            YLims = get(gca, 'YLim');
            
                        
            % plot lines at significant event tmes            
            for iEvent = 1:length(plot_obj.significant_event_times)
               line([plot_obj.significant_event_times(iEvent), plot_obj.significant_event_times(iEvent)], [YLims(1) YLims(2)], 'color', [0 0 0])
            end
            
            
            % plot shaded regions indicating where events occurred
            for iEvent = 1:length(plot_obj.significant_event_regions)
                
                start_time = [plot_obj.significant_event_regions{iEvent}(1)];
                end_time = [plot_obj.significant_event_regions{iEvent}(2)];
                YLims = get(gca, 'YLim');
                
                if isscalar(plot_obj.significant_event_region_alphas), 
                    curr_alpha = plot_obj.significant_event_region_alphas;
                elseif ismatrix(plot_obj.significant_event_region_alphas)
                    curr_alpha = plot_obj.significant_event_region_alphas(iEvent);
                else
                    error('significant_event_region_alphas must be either a scalar or a vector');
                end

                if ismatrix(plot_obj.significant_event_region_colors)                   
                    curr_color = plot_obj.significant_event_region_colors;
                elseif iscell(plot_obj.significant_event_region_colors)
                    curr_color = plot_obj.significant_event_region_colors{iEvent};
                else
                    error('significant_event_region_colors must be either a vector or a cell array of vectors');
                end
                
                 fill([start_time end_time end_time start_time],[YLims(1) YLims(1) YLims(2) YLims(2)], curr_color, 'EdgeColor', curr_color, 'EdgeAlpha', curr_alpha, 'FaceAlpha', curr_alpha);
            
            end
            
            
            if ~isempty(plot_obj.legend_names)
                legend(plot_obj.legend_names)
            end
            
            
            % make the tick marks point outward, only have axes on left and bottom, and set the paper position mode to auto 
            box off
            set(gca, 'tickdir', 'out')
            set(gcf, 'PaperPositionMode', 'Auto');  % will print the figure as it is displayed
            
            
        end

        
        
        
    end

end   % end classdef








