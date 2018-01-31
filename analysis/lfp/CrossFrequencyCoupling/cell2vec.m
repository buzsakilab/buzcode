function colvec = cell2vec(cel, inds_or_str, inds2)
% colvec = cell2vec(cel, inds)
%
% Created by EWS, September 2012

if (nargin < 2) || (isnumeric(inds_or_str) && ~isempty(inds_or_str) && (inds_or_str(1) < 1))
    inds = [];
elseif ischar(inds_or_str)
    if strcmpi(inds_or_str, 'all')
        inds = cell(length(cel),1);
        for i=1:length(cel)
            inds{i} = 1:length(cel{i});
        end
    elseif strcmpi(inds_or_str, 'ends')
        inds = cell(length(cel),1);
        for i=1:length(cel)
            inds{i} = [1 length(cel{i})];
        end
    elseif strcmpi(inds_or_str, 'except')
        inds = cell(length(cel),1);
        for i=1:length(cel)
            inds{i} = setdiff(1:length(cel{i}), inds2);
        end
    elseif strcmpi(inds_or_str, 'include')
        inds = cell(length(cel),1);
        for i=1:length(cel)
            inds{i} = intersect(1:length(cel{i}), inds2);
        end
    else
        error('invalid specifier string')
    end
elseif (length(inds_or_str) > 1) && (length(inds_or_str) < length(cel))
    error('inds argument must be either length 1 or same as cell argument')
elseif ~iscell(inds_or_str) && (length(inds_or_str) == 1)
    inds = ones(length(cel),1)*inds_or_str;
else
    inds = inds_or_str;
end

colcel = cel(:);
for i=1:length(colcel)
    colcel{i} = colcel{i}(:);
end
if isempty(inds)
    colvec = cell2mat(colcel);
elseif iscell(inds)
    colvec = [];
    for i=1:length(colcel)
        colvec = [colvec; colcel{i}(inds{i})];
    end
else
    colvec = zeros(length(colcel),1);
    for i=1:length(colcel)
        colvec(i) = colcel{i}(inds(i));
    end
end
