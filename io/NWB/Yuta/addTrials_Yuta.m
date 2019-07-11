function nwb = addTrials_Yuta(xml,nwb)
    %% THIS SECTION LOADS THE INFO FROM *RUN.MAT - YUTA - NON STANDARD FORMAT

    %  NOTE: Assumes the presence of a single *Run.mat and *RunInfo.mat file. If more are present, the code needs to change to store everything in: nwb.processing

    %% Add Trials

    % This file holds a matrix with the trial info
    trialsFile = dir([xml.folder_path filesep '*Run.mat']);

    if ~isempty(trialsFile)

        trialsInfo = load(fullfile(trialsFile.folder, trialsFile.name));
        the_field = fieldnames(trialsInfo);
        trialsInfo = trialsInfo.(the_field{1});

        colname = the_field{1};

        % This file holds a cell array with the labels for the matrix above (...)
        runFile = dir([xml.folder_path filesep '*RunInfo.mat']);
        runInfo = load(fullfile(runFile.folder,runFile.name));
        the_field = fieldnames(runInfo);
        runInfo = runInfo.(the_field{1});


        start_time_trials = types.core.VectorData('description','Starting timepoint of Each Trial','data',trialsInfo(:,1));
        stop_time_trials  = types.core.VectorData('description','Ending timepoint of Each Trial','data',trialsInfo(:,2));

        id_trials = types.core.ElementIdentifiers('data',int64(0:size(trialsInfo,1)-1)');

        conditions_trials = cell(size(trialsInfo,1),1);
        for iTrial = 1:size(trialsInfo,1)
            conditions_trials{iTrial} = runInfo{find(trialsInfo(iTrial,3:4))+2};
        end


        nwb.intervals_trials = types.core.TimeIntervals('start_time',start_time_trials,'stop_time',stop_time_trials,...
                                                   'colnames',colname,...
                                                   'description',['experimental trials from ' colname 'RunInfo.mat file'],'id',id_trials);

        nwb.intervals_trials.vectordata.set('both_visit', types.core.VectorData('description', 'Both visit condition', 'data', trialsInfo(:,7)));
        nwb.intervals_trials.vectordata.set('condition', types.core.VectorData('description', 'Condition Label', 'data', conditions_trials));
        nwb.intervals_trials.vectordata.set('error_run', types.core.VectorData('description', 'Error run', 'data', trialsInfo(:,5)));
        nwb.intervals_trials.vectordata.set('stim_run', types.core.VectorData('description', 'Stimulation run condition', 'data', trialsInfo(:,6)));

        disp('Yuta Trial info added..')

    else
        disp('No Yuta Trial info present in the selected directory..')
    end


end