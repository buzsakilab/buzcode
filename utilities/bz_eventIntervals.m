function eventInt = bz_eventIntervals(events,intervals)
%EVENTINT Separate given buzcode event struct by given time intervals
%   events - buzcode event structure containing pulses and 
%   intervals - set of time intervals
%   Requires: first field of events struct must be timestamps
%   Jacob Kerr, 2020

    loadStruct = load(events);
    load(intervals);
    
    loadFields = fieldnames(loadStruct);
    eventStruct = getfield(loadStruct,loadFields{1,1});
    fields = fieldnames(eventStruct);
    
    for i = 1:size(intervals)
        tmp = InIntervals(getfield(eventStruct,fields{1,1}),intervals(i,:)); % this function needs to be compiled
        timestampField = getfield(eventStruct,fields{1,1});
        eventInt{i}.timestamps = timestampField(tmp==1,:);
        for j = 2:size(fields)
            if size(getfield(eventStruct, fields{j,1}),1) == size(getfield(eventStruct,fields{1,1}),1)
                tempArr = getfield(eventStruct, fields{j,1});
                eventInt{i}.(fields{j,1}) = tempArr(tmp==1,:); 
            elseif size(getfield(eventStruct, fields{j,1}),2) == 2
                newTmp = InIntervals(getfield(eventStruct,fields{j,1}),intervals(i,:));
                tempArr = getfield(eventStruct, fields{j,1});
                eventInt{i}.(fields{j,1}) = tempArr(newTmp==1,:); 
            else 
                eventInt{i}.(fields{j,1}) = getfield(eventStruct, fields{j,1});
            end
        end
        %clear tmp1; 
    end
end
