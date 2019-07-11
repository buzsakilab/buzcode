function xml = GetXMLInfo(folder_path)

    %% This uses an .xml importer downloaded from MathWorks - File Exchange
    %  https://www.mathworks.com/matlabcentral/fileexchange/28518-xml2struct
    %  The loadxml from the Buzcode repo gave errors
    
    [previous_path, name] = fileparts(folder_path);

    all_files_in_folder = dir(folder_path);

    iXML = [];
    for iFile = 1:length(all_files_in_folder)
        if strfind(all_files_in_folder(iFile).name,'.xml')
            iXML = [iXML iFile];
        end
    end
    if isempty(iXML)
        error 'There are no .xml files in this folder'
    elseif length(iXML)>1
        error 'There is more than one .xml in this folder'
    end

    xml               = xml2struct([folder_path filesep all_files_in_folder(iXML).name]);
    xml               = xml.parameters;
    xml.folder_path   = folder_path;
    xml.name          = name;
end