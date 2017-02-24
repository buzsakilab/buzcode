function MClust_ApplyConvexHulls()

% MClust_ApplyConvexHulls
%
% Takes mcconvexhull clusters from a .clusters file and applies the hulls to the 
% current data set so as to create a new set of MClust_Clusters
%
% ADR 2008
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.5A.
% Version control M3.5A.


global MClust_Clusters
global MClust_TTdn MClust_TTfn

% Checks
	
if isempty(MClust_TTfn)
    errordlg('Select a TT and calculate features first before loading clusters.',...
        'MClust error', 'modal');
    return;
end

% Find Cluster file
[fn,dn] = uigetfile(fullfile(MClust_TTdn, [MClust_TTfn '.clusters']));
fn = fullfile(dn, fn);
if isempty(fn)
    errordlg('No clusters file found.', 'MClust error', 'modal');
    return;
end
temp = load(fn, '-mat');

% For each cluster, match, create, and append
nAdded = 0;
for iC = 1:length(temp.MClust_Clusters)
    C = temp.MClust_Clusters{iC};
    if isa(C, 'mcconvexhull')
        % convert
        newC = mcconvexhull(GetName(C), C, 'convert');
        if ~isempty(newC)
            % append to MClust_Clusters
            MClust_Clusters{end+1} = newC;
            MClust_Colors((length(MClust_Clusters)+1),:) = temp.MClust_Colors(iC+1,:);
            nAdded = nAdded+1;
        end
    end
end
        
% return
msgbox([num2str(nAdded) ' clusters appended to clusters list.'], 'MClust msg');

