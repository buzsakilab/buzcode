% FileExists(FileName)
%
% returns 1 if the file exists, 0 otherwise

function exists = FileExists_ss(FileName)

fp = fopen(FileName, 'r');

if (fp==-1)
	exists = 0;
else
	fclose(fp);
	exists = 1;
end;