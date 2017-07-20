function AlignTsp2Whl_All(varargin)

% USAGE
%     AlignTsp2Whl_All('parameters','values')
%     
% Inputs:
%    SIMPLEST THING TO DO: call AlignTsp2Whl_All()
%
%    WARNING: it is assumed that the root folder name and the final merge name are the same. If not, use
%    the 'mergename' option. Typically the root folder is something like
%    /xxx/yyy/AnimalZZZ-DayA/ and the merged data are in /xxx/yyy/AnimalZZZ-DayA/AnimalZZZ-DayA/
%
% - options:
%     'mergename': if the merged data are not in 'fbasename/fbasename' but in
%     'fbasename/mergename', then indicate with this option the name of the
%     merged data directory
%     'merge': if set to 1 (default 1) the function will merge the whl
%     dat into one single whl file 'fbasename/mergename/mergename.whl'.
%     'colorIx': a 2 value vector [colFront colRear] which defines which
%     color from [R G B] is the front and rear LEDs. Default: [1 3]*
%     'folderNames': a cell array of filenames in the order they will be
%     merged. If omitted, the program looks at the merged files in the xml
%     file contained in '/dir1/dir2/.../fbasename/mergename/mergename.xml'
%     'fromMpg': if 1 (default 0) will reconstruct the tracking from the
%     video
% 
% Dependencies: xml_toolbox (easy to find on the web)
% 
% Adrien Peyrache, 2012

merge=1;
colorIx=[1 3];
tspFiles = {};
tspFilesIx = [];
frommpg = 0;
[dummy fbasename dummy] = fileparts(pwd);

mergename = fbasename;

for i = 1:2:length(varargin),
  if ~isa(varargin{i},'char'),
    error(['Parameter ' num2str(i+3) ' is not a property (type ''help LoadBinary'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'mergename',
      mergename = varargin{i+1};
      if ~ischar(mergename)
        error('Incorrect value for property ''mergename''. Should be char array.');
      end 
    case 'merge',
      merge = varargin{i+1};
      if merge ~= 1 || merge ~= 0
        error('Incorrect value for property ''merge''. Should be logical.');
      end
    case 'colorix',
      colorIx = varargin{i+1};
      if length(colorIx)~=2
        error('Incorrect value for property ''colorIx''. Should be 3 element vector of logical values.');
      end
    case 'folderNames',
      tspFiles = varargin{i+1};
      isachar = 1;
      for fix=1:length(tspFiles),isachar=isachar*isa(tspFiles{fix},'char');end
      if ~isa(tspFiles,'cell') || ~isachar
        error('Incorrect value for property ''tspFiles''. Should be a cell array of filenames');
      end
     case 'frommpg',
         frommpg = varargin{i+1};
         if ~isnumeric(frommpg)
            error('Incorrect value for property ''frommpg''. Should be 0 pr 1');
         end
     case 'tspFilesix',
         tspFilesIx = varargin{i+1};
  end
end

% if isempty(tspFiles)
%     xmldata = xml_load([fbasename '.xml']);
%     segnames = xmldata.multiFileProcessing.files;
% 
%     if isempty(tspFilesIx)
%         tspFilesIx = (1:length(segnames));
%     else
%         tspFilesIx = tspFilesIx(:)';
%     end
%     tspFiles = cell(length(tspFilesIx),1);
%     for ii=1:length(tspFilesIx)
%         tspFiles{ii} = segnames(tspFilesIx(ii)).fileBaseName;
%     end
% end

tsptemp = dir('*.tsp');%get list of tsp files (to structure);
for a = 1:length(tsptemp)
    tspFiles{a}=tsptemp(a).name(1:end-4);
end

fprintf('Creating whl files...\n')
for fix=1:length(tspFiles)
    fname = tspFiles{fix};
    fprintf('\t%s\n',fname)
%     try
        AlignTsp2Whl(fname,colorIx,frommpg);    
%     catch
%         keyboard
%     end
end

if merge
    fprintf('Merging files...\n')
    whlall = [];
    for fix=1:length(tspFiles)
        fname = tspFiles{fix};
        whl = load([fname '.whl']);
        whlall = [whlall;whl];
    end
    dlmwrite([mergename '.whl'],whlall,'delimiter','\t');
end
