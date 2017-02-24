function C = List2Cell(filename)
% C = List2Cell(filename)
% INPUTS
% filename: a file containing a list of strings (one per line)
% 
% OUTPUTS
% C: a cell array containing the strings


fid = fopen(filename, 'r');
if(fid == -1)
  error('Could not open file');
end

S = fscanf(fid, '%c');

fclose(fid);

C = cell(0,1);
n = 1;

while(S)
  [d, S] = strtok(S);
  if ~isempty(d)
    if d(1) ~= '#'
      C{n} = d;
      n = n + 1;
    end
  end
end

if ~isempty(C)
  C = C(:);
else
  C = cell(0,1);
end
