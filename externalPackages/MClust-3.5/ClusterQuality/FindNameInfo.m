function NameInfo = FindNameInfo(filename)

% NameInfo = FindNameInfo(filename);
%
% From a filename (spike file, video tracker file, events file, etc), returns 
% an object whose fields contain the rat's identity, the session and 
% the tetrode and cluster number if available.
%
% INPUTS: filename: expects a name including a path, unless the name is in the format RXXX-YYYY-MM-DD-TTtt-cc.*
%                   expects rat ID in RXXX, session date in YYYY-MM-DD
%
% OUTPUTS: NameInfo.RatID = rat's identification number
%          NameInfo.FileName = original filename passed in
%          NameInfo.SessionID = RXXX-YYYY-MM-DD
%          NameInfo.Tetrode = number of the tetrode
%          NameInfo.Cluster = number of the cluster
%          NameInfo.Location = directory of filename
%
% ncst 7 June 02
% status PROMOTED
%
% modified ncst 21 Jun 02 to handle .t files starting with Sc (new Neuralynx naming format)

NameInfo = [];
FileName = filename;
RatID = [];
SessionID = [];
TetrodeID = [];
Tetrode = [];
Cluster = [];
Location = [];
TetrodePrefix = [];

% Initialize output variable 
NameInfo.FileName = filename;
NameInfo.Location = [Location filesep];
NameInfo.RatID = RatID;
NameInfo.SessionID = SessionID;
NameInfo.TetrodePrefix = TetrodePrefix;
NameInfo.Tetrode = Tetrode;
NameInfo.Cluster = Cluster;
NameInfo.TetrodeID = TetrodeID;
NameInfo.CellID = TetrodeID;

if isa(filename,'cell')
    old_filename = filename;
    clear filename
    filename = old_filename{1};
    clear old_filename
end

[p n e] = fileparts(filename);

if isempty(p)
    old_fn = filename;
    filename = findfiles(filename);
    if isempty(filename)
        disp([old_fn ' does not exist'])
        return
    end
    old_filename = filename;
    clear filename
    filename = old_filename{1};
    clear old_filename
end

[p n e] = fileparts(filename);

Location = p;

Rs = findstr(p,'R');
if ~isempty(Rs)
    RatID = p(Rs(end):Rs(end) + 3);
    SessionID = p(Rs(end):Rs(end) + 14);
else
    Rs = findstr(n,'R');
    if ~isempty(Rs)
        RatID = n(Rs(end):Rs(end) + 3);
        SessionID = n(Rs(end):Rs(end) + 14);
    end
end

TTs = [findstr(n,'TT') findstr(n,'Sc')];
UndScr = findstr(n,'_');
if ~isempty(TTs)
    if ~isempty(UndScr)
        Tetrode = n(TTs(end)+2:UndScr-1);
        Cluster = n(UndScr+1:end);
        if ~isempty(findstr(n,'Sc'))
            TetrodePrefix = 'Sc';
        else
            TetrodePrefix = 'TT';
        end
    else
        Tetrode = n(TTs(end) + 2:TTs(end)+3);
        Cluster = n(TTs(end) + 5:end);
        if exist([Location filesep 'Sc' num2str(str2num(Tetrode)) '.ntt'])
            TetrodePrefix = 'Sc';
        else
            TetrodePrefix = 'TT';
        end
    end
    TetrodeID = [SessionID '-TT' Tetrode '-' Cluster];
end

% if ~isempty(Cluster)
%     if str2num(Cluster) < 10
%         Cluster = ['0' Cluster];
%     end
% end

NameInfo.FileName = filename;
NameInfo.Location = [Location filesep];
NameInfo.RatID = RatID;
NameInfo.SessionID = SessionID;
NameInfo.TetrodePrefix = TetrodePrefix;
NameInfo.Tetrode = Tetrode;
NameInfo.Cluster = Cluster;
NameInfo.TetrodeID = TetrodeID;
NameInfo.CellID = TetrodeID;
NameInfo.CellID = TetrodeID;
    
