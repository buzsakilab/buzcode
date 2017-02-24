function S = GetListOfFiles(tstring)
%
% S =  GetListOfFiles(tstring)
%
% inp: a string to be used as a regexp for identifying files within a
% directory p. es. ('trail')
% 
% out: a cell array of strings containing the file names
% 
% batta 1999
% version: ignota
% status: BASTACHECAMMINA


if ~isa(tstring, 'char')
  error('GetListOfTfiles: argument must be a string');
end

if length(tstring) < 1
  error('GetListOfTfiles: argument must be non empty');
end
%regexp = [ '*' tstring '*'];
regexp = tstring;
lsstring = ls(regexp);

if length(lsstring) < 1
  return;
end

while ~isempty(findstr(lsstring, 'ls:'))
  errmsgstart = findstr(lsstring, 'ls:');
  errmsgstart = errmsgstart(1);
  errmsgstop = findstr(lsstring, 'directory');
  errmsgstop = errmsgstop(1) +   8;
  lsstring = lsstring([1:errmsgstart-1 errmsgstop+2:length(lsstring)]);
end



n_elem = 1;
start = 1;
stop = 1;
S = cell(1,0);
while stop <= length(lsstring)
  while ~isspace(lsstring(stop)) % & stop <= length(lsstring)
    stop = stop + 1;
  end
  S{n_elem} = lsstring(start:stop-1);
  n_elem = n_elem + 1;
  start = stop + 1;
  if start > length(lsstring)
    break
  end
  
  while(isspace(lsstring(start)))
    start = start + 1;
    end

  stop = start;
end
