function sizeTree = folderSizeTree(directory)
%FOLDERSIZETREE - Calculate total size of sub-folders recursively
%
%	usage:	sizeTree = folderSizeTree(directory)
%	The function scans all subfolders of given "directory", and returns
%	a list of all subfolders and their size.
%	The size of a folder is the sum of all its files size + sum of all its
%	subfolders size.
%
%	returned value is a struct of 3 coordinated cells:
%		sizeTree.name => list of all subfolders.
%		sizeTree.size => size of each subfolder.
%		sizeTree.level => level of each subfolder. 0 is given folder.
%								1 is one level deep etc.
%
%	Copyright: Yanai Ankri, November 2010



sizeTree.name = {};
sizeTree.size = {};
sizeTree.level = {};
size = 0;

%scan folder
list = dir(directory);
for i=1:length(list)
	fileName = list(i).name;
	
	if isequal(fileName, '.') || isequal(fileName, '..'),		continue;		end

	%if "file" is subfolder, scan it recursively
	if list(i).isdir
		fullName = fullfile(directory, fileName);
		sizeSub = folderSizeTree(fullName);

		%add the result to the main list, and increase the level value
		sizeI = 0;
		for j=1:length(sizeSub.name)
			sizeTree.name{end+1} = sizeSub.name{j};
			sizeTree.size{end+1} = sizeSub.size{j};
			sizeTree.level{end+1} = sizeSub.level{j}+1;
			if sizeSub.level{j} == 0	%to prevent summing subfolders
				sizeI = sizeI + sizeSub.size{j};
			end
		end
	else
		%"file" is realy a file (not subfolder)
		sizeI = list(i).bytes;
	end

	%total size of this folder
	size = size + sizeI;
end

sizeTree.name{end+1} = directory;
sizeTree.size{end+1} = size;
sizeTree.level{end+1} = 0;


end

		
