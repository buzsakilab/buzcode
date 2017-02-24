function S = LoadSpikes(tfilelist)

% S = LoadSpikes(tfilelist)
% inp: tfilelist is a cellarray of strings, each of which is a
% 	tfile to open.  Note: this is incompatible with version unix3.1.
% out: Returns a cell array such that each cell contains a ts 
% 	object (timestamps which correspond to times at which the cell fired)

% ADR 1998
%  version L4.0
%  status: PROMOTED

% Francesco P. Battaglia 2004
% version modified to use the tsdArray class  
  
  
%-------------------
% Check input type
%-------------------
if ~isa(tfilelist, 'cell')
   error('LoadSpikes: tfilelist should be a cell array.');
end

nFiles = length(tfilelist);

%--------------------
% Read files
%--------------------

%fprintf(2, 'Reading %d files.', nFiles);

% for each tfile
% first read the header, the read a tfile 
% note: uses the bigendian modifier to ensure correct read format.

S = cell(nFiles, 1);
for iF = 1:nFiles
%   DisplayProgress(iF, nFiles, 'Title', 'LoadSpikes');
  tfn = tfilelist{iF};
  if ~isempty(tfn)
    tfp = fopen(tfn, 'rb','b');
    if (tfp == -1)
      warning([ 'Could not open tfile ' tfn]);
    end
    
    skipHeader(tfp);    
    S{iF} = fread(tfp,inf,'uint32');	%read as 32 bit ints
    S{iF} = ts(S{iF}, 'fixOrder', 1);
    setName(S{iF}, tfn);
    fclose(tfp);
 
  end 		% if tfn valid
end		% for all files
%fprintf(2,'\n');

S = tsdArray(S);


function skipHeader(fp)

done = false;
n = length('%%ENDHEADER');
while ~done 
    s = fgets(fp);
    if strncmp(s, '%%ENDHEADER', n)
        done = true;
    end
end



