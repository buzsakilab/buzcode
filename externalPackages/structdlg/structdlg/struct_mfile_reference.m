function S = struct_mfile_reference(mfile,struct_name,strs)
% S = struct_mfile_reference(mfile,struct_name,strs)

% AF 11/1/01

pattern = [struct_name '.'];
if (~isempty(mfile))
   strs = grep(mfile,pattern);
end
S = [];

for i = 1:length(strs)
   inds = strfind(strs{i},pattern);
   for j = inds(:)'
      start_field = j+length(pattern);
      k = start_field+1;
      while ((k <= length(strs{i})) & (isvarname(strs{i}(start_field:k))))
         k = k+1;
      end
      field_name = strs{i}(start_field:k-1);
      if (~isempty(field_name))
         eval(['S.' field_name ' = struct_mfile_reference(''' mfile ''',''' struct_name '.' field_name ''',strs);']);
      end
   end
end
