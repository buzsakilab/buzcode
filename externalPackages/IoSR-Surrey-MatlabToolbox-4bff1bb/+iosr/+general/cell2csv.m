function cell2csv(C,filename)
%CELL2CSV Output a cell array to a CSV file
% 
%   IOSR.GENERAL.CELL2CSV(C,FILENAME) writes the cell array C to a CSV file
%   specified by filename.

%   Copyright 2016 University of Surrey.

    %% make sure filename has a .csv extension

    [pathstr, name, ext] = fileparts(filename);
    if ~strcmp(ext,'.csv')
        filename = [pathstr filesep name '.csv'];
    end

    %% write data

    [nrows,ncols]= size(C);

    % trap and remove cell arrays of strings
    IX = cell2mat(cellfun(@iscell,C,'UniformOutput',false));
    C(IX) = cellfun(@char,C(IX),'UniformOutput',false);

    % create file
    fid = fopen(filename, 'w');

    % create format string
    format = cell(1,ncols);
    ind = cellfun(@ischar,C(2,:));
    format(ind) = {'%s,'};
    format(~ind) = {'%f,'};
    format = wrap(format);

    if any(cellfun(@ischar,C(1,:))~=cellfun(@ischar,C(2,:)))
        % the cell array has a heading line
        header = repmat({'%s,'},[1,ncols]);
        header = wrap(header);
        start = 2;
        fprintf(fid, header, C{1,:});
    else
        % the cell array does not have a heading line
        start = 1;
    end

    for row=start:nrows
        fprintf(fid, format, C{row,:});
    end

    fclose(fid);

end

function y = wrap(x)
%WRAP Concatenate cells and add end EOL

    y = cell2mat(x);
    i = strfind(y,',');
    y = [y(1:i(end)-1) '\n'];

end
