function compilefma

% compilefma - Compile C/C++ functions in FMAToolbox

% Move to the 'FMAToolbox' directory, compile, and move back
path = fileparts(which('FMAToolbox'));
currentDir = pwd;
cd(path);
processDir;
cd(currentDir);

function processDir

files = dir;
for i = 1:length(files),
	if strcmp(files(i).name,'.') || strcmp(files(i).name,'..'), continue; end
	if files(i).isdir,
		% Recursively process subdirectories
		cd(files(i).name);
		processDir;
		cd('..');
	end
	% Compile C/C++ files
	[unused,name,ext] = fileparts(files(i).name);
	if strcmp(lower(ext),'.c') || strcmp(lower(ext),'.cpp'),
		mex(files(i).name);
	end
end
