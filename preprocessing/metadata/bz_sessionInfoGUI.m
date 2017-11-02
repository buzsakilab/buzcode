function [ sessionInfo ] = bz_sessionInfoGUI( sessionInfo )
%A gui for editing buzcode sessionInfo structures.  Right now it just adds
%region data, but will be expanded soon.
%%

%% add region window
%Get the spike groups into selection format
for ss = 1:sessionInfo.spikeGroups.nGroups
    spikegroups.(['sg',num2str(ss)]) = { {'{0}','1'},...
        ['SpikeGroup',num2str(ss),...
        ' (',num2str(sessionInfo.spikeGroups.groups{ss}(1)),', ... ,',...
        num2str(sessionInfo.spikeGroups.groups{ss}(end)),')']};
end

%Run the user input
newstruct.anotherregion = true; numregions = 0;
while newstruct.anotherregion
    %Make the template structure
    tempstruct.loadexisting = { {'uigetfile(''*.sessionInfo.mat'')'} };
    tempstruct.(['regionname',num2str(numregions)]) = ...
        { 'HPC/CTX/other?' ,['Region ',num2str(numregions),' Tag']};
    tempstruct.(['regionchans',num2str(numregions)]) = ...
        {[] ,'Channels (0-idx)'};
    tempstruct.Select_SpikeGroups = spikegroups;
    tempstruct.anotherregion = { {'{0}','1'} 'Add Another Region'};
    
    %Run StructDlg to get user input
    newstruct = StructDlg(tempstruct,'Add Region Info.');
    
    %For loading from an existing sessionInfo
    if newstruct.loadexisting
        oldsessionInfo = load(newstruct.loadexisting);
        sessionInfo.region = oldsessionInfo.sessionInfo.region;
        return %return the sessionInfo with regions from old file...
    end
        
    %Extract spike groups
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
end

%Add the regions to the sessionInfo structure
sessionInfo.region = cell(1,sessionInfo.nChannels);
for rr = 1:numregions
    regionname = newstruct.(['regionname',num2str(rr-1)]);
    regionchans = newstruct.(['regionchans',num2str(rr-1)]);
    sessionInfo.region(ismember(sessionInfo.channels,regionchans)) = {regionname};
end

end

