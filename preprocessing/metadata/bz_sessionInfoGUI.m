function [ sessionInfo ] = bz_sessionInfoGUI(  sessionInfo,MODE )
%A gui for editing buzcode sessionInfo structures.  Right now it just adds
%region data, but will be expanded soon.
%
%INPUTS
%   sessionInfo     a buzcode sessionInfo structure (see wiki) to edit
%   MODE            (optional) - a field to edit
%                   options: 'Regions','Bad Channels'
%OUTPUTS
%   sessionInfo     the edited sesionInfo structure
%
%DLevenstein 2017
%NOTE: this GUI is v much in prgress (beta)... please improve!
%%
if ~exist('MODE','var')
    MODE = 'none';
end


fieldoptions = {'Regions','Bad Channels','quit'};

%Get the spike groups into selection format
for ss = 1:sessionInfo.spikeGroups.nGroups
    spikegroups.(['sg',num2str(ss)]) = { {'{0}','1'},...
        ['SpikeGroup',num2str(ss),...
        ' (',num2str(sessionInfo.spikeGroups.groups{ss}(1)),', ... ,',...
        num2str(sessionInfo.spikeGroups.groups{ss}(end)),')']};
end

%% Loop For user input modes
while ~strcmp(MODE,'quit')
switch MODE
    
    case 'none'     %%Modeselektor
        [s,v] = listdlg('PromptString','Edit field:',...
                'SelectionMode','single',...
                'ListString',fieldoptions);
            if v==0
                MODE = 'quit';
            else
                MODE = fieldoptions{s};
            end

            
    %%MODE: add region
    case 'Regions'  

        
        clear tempstruct
        tempstruct.Load_From_Existing = { {'uigetfile(''*.sessionInfo.mat'')'} };
        
        %Do some regions already exist?
        if isfield(sessionInfo,'region')
            existingregions = unique(sessionInfo.region);
            existingregions(strcmp(existingregions,''))=[];
            numregions = length(unique(existingregions));
            for rr = 1:numregions
                tempstruct.(['regionname',num2str(rr-1)]) = ...
                    existingregions{rr};
                tempstruct.(['regionchans',num2str(rr-1)]) = ...
                    sessionInfo.channels(strcmp(sessionInfo.region,existingregions{rr}));
           % existingbadchans = sessionInfo.region;
            end
        else
            numregions = 0;
        end
        


        %Run the user input
        newstruct.anotherregion = true; 
        while newstruct.anotherregion
            %Make the template structure
            tempstruct.(['regionname',num2str(numregions)]) = ...
                { 'HPC/CTX/other?' ,['Region ',num2str(numregions),' Tag']};
            tempstruct.(['regionchans',num2str(numregions)]) = ...
                {[] ,'Channels (0-idx)'};
            tempstruct.Select_SpikeGroups = spikegroups;
            tempstruct.anotherregion = { {'{0}','1'} 'Add Another Region'};

            %Run StructDlg to get user input
            newstruct = StructDlg(tempstruct,'Add Region Info.');
            
            %For Cancel
            if isempty(newstruct)
                break 
            end
            
            %Loading from an existing sessionInfo
            if ~strcmp(newstruct.Load_From_Existing,'*.sessionInfo.mat')
                oldsessionInfo = load(newstruct.Load_From_Existing);
                sessionInfo.region = oldsessionInfo.sessionInfo.region;
                return %return the sessionInfo with regions from old file...
            end

            %Extract spike groups from user input
            for ss = 1:sessionInfo.spikeGroups.nGroups
                if newstruct.Select_SpikeGroups.(['sg',num2str(ss)])
                    newstruct.(['regionchans',num2str(numregions)]) = ...
                        [newstruct.(['regionchans',num2str(numregions)]),...
                        sessionInfo.spikeGroups.groups{ss}];
                end
            end

            %Reset the template structure and count up the regions
            tempstruct = rmfield(newstruct,{'anotherregion','Select_SpikeGroups'});
            numregions = numregions+1;
        end     %End the user input loop

        %Add the regions to the sessionInfo structure
        if ~isempty(newstruct) %For user:cancel
            sessionInfo.region = repmat({''},1,sessionInfo.nChannels);
            for rr = 1:numregions
                regionname = newstruct.(['regionname',num2str(rr-1)]);
                regionchans = newstruct.(['regionchans',num2str(rr-1)]);
                sessionInfo.region(ismember(sessionInfo.channels,regionchans)) = {regionname};
            end
        end
        MODE = 'none'; %Return to selection
        
    %MODE: Bad Channels    
    case 'Bad Channels'
        %Do some badchannels already exist?
        if isfield(sessionInfo,'badchannels')
            existingbadchans = sessionInfo.badchannels;
        else
            existingbadchans = [];
        end
        
        %Run the user input
        clear tempstruct
        %Make the template structure
        tempstruct.Load_From_Existing = { {'uigetfile(''*.sessionInfo.mat'')'} };
        tempstruct.badchans = ...
            {existingbadchans ,'Channels (0-idx)'};
        tempstruct.Select_SpikeGroups = spikegroups;
        newstruct = StructDlg(tempstruct,'Add Bad Channels');
        
        if ~isempty(newstruct) %For user: cancel
            %Loading from an existing sessionInfo
            if ~strcmp(newstruct.Load_From_Existing,'*.sessionInfo.mat')
                oldsessionInfo = load(newstruct.Load_From_Existing);
                sessionInfo.badchannels = oldsessionInfo.sessionInfo.badchannels;
                return %return the sessionInfo with badchannels from old file...
                %       NOTE: should not return, but should go back to
                %       ModeSelektor?
            end

            %Extract spike groups from user input
            for ss = 1:sessionInfo.spikeGroups.nGroups
                if newstruct.Select_SpikeGroups.(['sg',num2str(ss)])
                    newstruct.badchans = ...
                        [newstruct.badchans,...
                        sessionInfo.spikeGroups.groups{ss}];
                end
            end

            %Add the badchannels to the sessionInfo structure
            sessionInfo.badchannels = newstruct.badchans;
        end
        MODE = 'none'; %Return to selection
        
    case 'quit'
        return
end
end

end


