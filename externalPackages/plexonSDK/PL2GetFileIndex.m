function  [pl2] = PL2GetFileIndex(filename)
% PL2GetFileIndex(filename): gets pl2 file index. File indexes are cached.
%        If cached file index for the specified filename exists, returns cached file index,
%       otherwise, calls PL2ReadFileIndex(filename)
%
% [pl2] = PL2GetFileIndex(filename)
%
% INPUT:
%   filename - pl2 file path. If empty string, will clear all cached file indexes. 
%
% OUTPUT:
%   pl2 - pl2 file index object (described in detail in PL2ReadFileIndex.m)
%       keeps pl2 file indexes cached in memory.
%       To clear indexes, use this command:
%       PL2GetFileIndex('')
%       or
%       clear PL2GetFileIndex
%


% this array of pl2 indexes is persisted between calls to this function
persistent PL2Indexes;

if nargin ~= 1
    error 'expected 1 input argument';
end

if isempty(filename)
    pl2 = 0;
    PL2Indexes={};
    return
end

if exist(filename,'file') ~= 2
    error('file does not exist');
end

if isempty(PL2Indexes)
    PL2Indexes={};
    PL2Indexes{1}.Index = PL2ReadFileIndex(filename);
    PL2Indexes{1}.Name = filename;
    pl2 = PL2Indexes{1}.Index;
    return;
else
    n = size(PL2Indexes,2);
    for i=1:n
        if strcmp(PL2Indexes{i}.Name, filename) == 1
            % we found Index with specified file Name
            if isfield(PL2Indexes{i}, 'Index')
                % if the index structure has pl2 index assigned, return it
                pl2 = PL2Indexes{i}.Index;
                return;
            else
                % read pl2 index from file, assign to persistent indexes and return
                PL2Indexes{i}.Index = PL2ReadFileIndex(filename);
                pl2 = PL2Indexes{i}.Index;
                return;
            end
        end
    end
    % we did not find persisted Index. 
    % read the index from file
    PL2Indexes{n+1}.Index = PL2ReadFileIndex(filename);
    PL2Indexes{n+1}.Name = filename;
    pl2 = PL2Indexes{n+1}.Index;
end

end

