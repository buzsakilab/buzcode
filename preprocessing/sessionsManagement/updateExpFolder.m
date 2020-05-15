
function updateExpFolder(inputFolder, outputFolder)
% Update experiment folder from Recording computer to Analysis computer
%
% USAGE 
%   updateExpFolder(inputFolder, outputFoder)
%
% INPUT
%   inputFolder     Experiment folder in recording computer, can be
%                       multiple folders
%   outputFolder    Experiment folder in analysis computer. Only one
%                       folder...
%
% Manu Valero-BuzsakiLab 2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
or = pwd;
% get session list codes from output folder
cd(outputFolder);
expName = strsplit(pwd,filesep);
expName = expName{end};
allRecOut = dir(strcat(expName,'_*'));
sessOutput(1,1:2) = NaN; 
for ii = 1:size(allRecOut,1)
    tmp = strsplit(allRecOut(ii).name,'_');
    sessOutput(ii,1) = str2num(tmp{2});
    if ~isempty(str2num(tmp{3}))
        sessOutput(ii,2) = str2num(tmp{3});
    else
        sessOutput(ii,2) = NaN;
    end
end

% if only one output folder is provided convert to cell
if ~iscell(inputFolder)
    ifol = inputFolder; clear inputFolder
    inputFolder{1} = ifol;
end

for jj = 1:length(inputFolder)
    fprintf(' Input folder %i of %i ...\n',jj,length(inputFolder));
    
    cd(inputFolder{jj});
    expNameInput = strsplit(pwd,filesep);
    expNameInput = expNameInput{end};
    if ~strcmp(expName,expNameInput)
        error('Experimental name does not match!!');
    end
    
    allRecInp = dir(strcat(expName,'_*'));
    for ii = 1:size(allRecInp,1)
        tmp = strsplit(allRecInp(ii).name,'_');
        recInput(ii,1) = str2num(tmp{2});
        if ~isempty(str2num(tmp{3}))
            recInput(ii,2) = str2num(tmp{3});
        else
            recInput(ii,2) = NaN;
        end
    end
    
    % folder not assigned to a session
    newExp = find(~(ismember(recInput(:,1), sessOutput(isnan(sessOutput(:,2)),1)))); 
    
    % discard already transfered folder
    newExp(ismember(recInput(newExp,2),sessOutput(:,2))) = [];
    % discard folders that are growing (active session)
    for ii = 1: length(newExp)
        cd(strcat(allRecInp(newExp(ii)).folder,filesep,allRecInp(newExp(ii)).name));
        f = dir('amplifier*.dat');
        f1 = fopen(f.name);
        fseek(f1, 0, 'eof');
        filesize(1) = ftell(f1);
        pause(0.1);
        fseek(f1, 0, 'eof');
        filesize(2) = ftell(f1);
        if (filesize(1) ~= filesize(2))
            newExp(ii) = 0;
        end
    end
    newExp(newExp==0) = [];
    cd(inputFolder{jj});

    for ii = 1:length(newExp)
        fprintf(' ** Copying %i file of %i ...\n',ii,length(newExp));
        copyfile(strcat(allRecInp(newExp(ii)).folder,filesep,allRecInp(newExp(ii)).name),...
            strcat(outputFolder,filesep,allRecInp(newExp(ii)).name));
    end
    clear expNameInput recInput allRecInp
end

cd(or);

end
