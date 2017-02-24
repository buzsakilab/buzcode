function [] = PL2Print(cellArrayOfStructures)
% PL2Print(cellArrayOfStructures): prints cell array of structures
%
% PL2Print(cellArrayOfStructures)
%
% INPUT:
%   cellArrayOfStructures - cell array of identical structures
%
% OUTPUT:
%   prints specified cell array

if nargin ~= 1
    error 'expected 1 input argument';
end

if numel(cellArrayOfStructures) == 0
    return
end 

names = fieldnames(cellArrayOfStructures{1});

% print field names
fprintf('[cell#]');
for i=1:numel(names)
    fprintf(' %18.18s', names{i});
end
fprintf('\n');

% print field values for each structure
for i=1:numel(cellArrayOfStructures)
    fprintf('%7d', i);
    theStructure = cellArrayOfStructures{i};
    for j=1:numel(names)
        x = theStructure.(names{j});
        switch class(x)
            case 'char'
                fprintf(' %18.18s', x);
            case 'double'
                if numel(x) == 0
                    fprintf('                   ');
                elseif numel(x) == 1
                    fprintf(' %18.10g', x);
                elseif numel(x) > 1 
                    fprintf(' %15.9g,..', x(1));
                end            
            otherwise
                    fprintf('                   ');
        end
    end
    fprintf('\n');
end
end
