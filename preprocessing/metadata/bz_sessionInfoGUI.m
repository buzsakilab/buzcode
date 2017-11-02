function [ sessionInfo ] = bz_sessionInfoGUI( sessionInfo )
%A gui for editing buzcode sessionInfo structures.  Right now it just adds
%region data, but will be expanded soon.
%%

%% add region window

clear tempstruct newstruct
newstruct.anotherregion = true;
numregions = 0;
while newstruct.anotherregion
    tempstruct.(['regionname',num2str(numregions)]) = { 'HPC/CTX/other?' ,...
        ['Region ',num2str(numregions),' Tag']};
    tempstruct.(['regionchans',num2str(numregions)]) = {[0:5] ,'Channels (0-idx)'};
    tempstruct.anotherregion = { {'{0}','1'} 'Add Another Region'};

    newstruct = StructDlg(tempstruct,'Add Region Info.');
    tempstruct = rmfield(newstruct,'anotherregion');
    numregions = numregions+1;
end

%Add the regions to the sessionInfo
parameters.region = cell(1,parameters.nChannels);
for rr = 1:numregions
    regionname = newstruct.(['regionname',num2str(rr-1)]);
    regionchans = newstruct.(['regionchans',num2str(rr-1)]);
    sessionInfo.region(ismember(sessionInfo.channels,regionchans)) = {regionname};
end
end

