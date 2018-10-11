function [str,cur_line] = struct2str(s,units,max_width,max_struct_elem, ident_str,total_ident_len,str,cur_line)
%

% AF 11/6/01

if (exist('units','var') ~= 1)
   units = struct([]);
end
if (exist('max_width','var') ~= 1)
   max_width  = Inf;
end
if (exist('max_struct_elem','var') ~= 1)
   max_struct_elem = 3;
end
if (exist('ident_str','var') ~= 1)
   ident_str  = '|';
end
if (exist('total_ident_len','var') ~= 1)
   total_ident_len  = 2;
end
if (exist('str','var') ~= 1)
   str = repmat({''},400,1);
   first_call = 1;
else
   first_call = 0;
end
if (exist('cur_line','var') ~= 1)
   cur_line = 0;
end
spacing = 2;

fnames = fieldnames(s);
fnames_lbl = build_labels(fnames,units);
max_lbl_width = size(char(fnames_lbl),2);
for i = 1:length(fnames)
   for j = 1:spacing-1
      cur_line = cur_line+1;
      str{cur_line} = ident_str;
   end
   cur_line = cur_line+1;
   str{cur_line} = ident_str;
   leading_spaces = repmat('-', 1, total_ident_len -length(ident_str)+max_lbl_width -length(fnames_lbl{i}));
   str{cur_line} = sprintf('%s%s%s: ', str{cur_line} ,leading_spaces, fnames_lbl{i});
   x = getfield(s,fnames{i});
   %% recursive call for sub-structures
   if (isstruct(x))
      new_ident_len = total_ident_len + max_lbl_width+2;
      new_ident_str = [ident_str repmat(' ',1,new_ident_len-2 - length(ident_str) - ceil(length(fnames_lbl{i})/2)) '|'];
      for xi = 1:min(length(x),max_struct_elem)
         if (isfield(units,fnames{i}))
            sub_units = getfield(units,fnames{i});
         else
            sub_units = struct([]);
         end
         [str,cur_line] = struct2str(x(xi),sub_units,max_width,max_struct_elem, new_ident_str,new_ident_len,str,cur_line);
         cur_line = cur_line+1;
         str{cur_line} = ident_str;
      end
      if (length(x) > max_struct_elem)
         dotted_str = [ident_str repmat(' ',1,new_ident_len-2 - length(ident_str) - ceil(length(fnames_lbl{i})/2)) ':'];
         for dot_i = 1:2
            cur_line = cur_line+1;
            str{cur_line} = dotted_str;
         end
      end
   else
      xstr = element2str(x,max_width-max_lbl_width-2);
      str{cur_line} = sprintf('%s%s', str{cur_line}, xstr);
   end
end
if (first_call)
   str = str(1:cur_line);
end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
%%                        SUB FUNCTIONS                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
function fnames_lbl = build_labels(fnames,units);
%
fnames_lbl = strrep(fnames,'_',' ');
f_units = fieldnames(units);
v_units = struct2cell(units);
for i = 1:length(f_units)
   if (ischar(v_units{i}) & ~isempty(v_units{i}))
      index = strmatch(f_units{i},fnames,'exact');
      if (~isempty(index))
         fnames_lbl{index} = strrep(v_units{i},'*',fnames_lbl{index});
         % fnames_lbl{index} = [fnames_lbl{index} ' (' v_units{i} ')'];
      end
   end
end
return;

% function fnames_lbl = build_labels(fnames,f_units);
% %
% fnames_lbl = strrep(fnames,'_',' ');
% for i = 1:min(length(fnames_lbl),length(f_units))
%    if (ischar(f_units{i}) & ~isempty(f_units{i}))
%       fnames_lbl{i} = [fnames_lbl{i} ' (' f_units{i} ')'];
%    end
% end
% return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xstr = element2str(x,max_width)
if (exist('max_width','var') ~= 1)
   max_width = Inf;
end
switch (class(x))
case 'char'
   if (length(x) < max_width-2)
      xstr = ['''' x ''''];
   else
      xstr = ['''' x(1:max_width-4) '...'''];
   end
   
case {'double' 'sparse'}
   if (isempty(x))
      xstr = '[]';
   elseif (ndims(x) > 2 | min(size(x)) >1 | length(x) >150)
      dims = size(x);
      xstr = ['[' sprintf('%dx',dims(1:end-1)) sprintf('%d',dims(end)) ' ' class(x) ']'];
   else 
      % x is a vector
      if (size(x,2) == 1)
         sep = ' ; ';
      else
         sep = ' ';
      end
      if (length(x) == 1)
         xstr = num2str(x);
      else
         xstr = ['[' num2str(x(1))];
         for ix = 2:length(x)
            xstr = [xstr sep num2str(x(ix))];
         end
         xstr = [xstr ']'];
      end
      if (length(xstr) > max_width)
         xstr = [xstr(1:max_width-4) '...]'];
      end
   end
   
case 'cell'
   xstr = '{';
   if (isempty(x))
      xstr = '{}';
   elseif (ndims(x) > 2 | min(size(x)) >1)
      dims = size(x);
      xstr = ['{' sprintf('%dx',dims(1:end-1)) sprintf('%d',dims(end)) ' cell}'];
   else 
      % x is a cell vector
      if (size(x,2) == 1)
         sep = ' ; ';
      else
         sep = ' ';
      end
      xstr = ['{' element2str(x{1},max_width/3)];
      for ix = 2:length(x)
         xstr = [xstr sep element2str(x{ix},max_width/3)];
      end
      xstr = [xstr '}'];
      if (length(xstr) > max_width)
         xstr = [xstr(1:max_width-4) '...}'];
      end
   end
   
case 'function_handle'
   xstr = element2str(['@' func2str(x)],max_width);
   
end
