function [ newstruct1, newstruct2 ] = bz_Matchfields( struct1,struct2,mode )
%[ newstruct1, newstruct2 ] = bz_Matchfields( struct1,struct2,mode )
%Takes two structures that may have different fields and gives them the
%same fields. Mode can either be 'remove' or 'add' (mismatching fields)
%
%If struct1 is a cell array of structues, and struct2 is [], then will
%return a structure array with matched fields
%%

%There's a better way to do this...
%allfields = cellfun(@(X) fieldnames(X),structs,'uniformoutput',false);
if iscell(struct1) & isempty(struct2)
    newstruct1 = struct1{1};
   for ss = 2:length(struct1)
       [newstruct1,newstruct2] = bz_Matchfields(newstruct1,struct1{ss},mode);
       newstruct1(ss) = newstruct2;
   end
   return
end


fieldnames1 = fieldnames(struct1);
fieldnames2 = fieldnames(struct2);

[commonfields,IdxFn1,IdxFn2] = intersect(fieldnames1,fieldnames2,'stable');

switch mode
    case 'remove'
        exclude1 = true(size(fieldnames1));exclude1(IdxFn1)=false;
        exclude2 = true(size(fieldnames2));exclude2(IdxFn2)=false;
        newstruct1 = rmfield(struct1,fieldnames1(exclude1));
        newstruct2 = rmfield(struct2,fieldnames2(exclude2));
        newstruct2 = orderfields(newstruct2,newstruct1);
        
    case 'add'
        
end
%%
% matchedstruct = structs{1};
% FIELDMISMATCH=false;
% for ss = 2:length(structs)
%     nextfields = fieldnames(structs{ss});
%     oldfields = fieldnames(matchedstruct);
%     newfields = setdiff(nextfields,oldfields);
%     
%     if ~isempty(newfields)
%         for ff = 1:length(newfields)
%             matchedstruct(1).(newfields{ff}) = []; 
%         end
%         FIELDMISMATCH=true;
%     end
% 
%     matchedstruct = orderfields(matchedstruct,structs{ss});   
%     
%     matchedstruct(ss) = structs{ss};
% end
%%
% cat(1,structin(:).(currentfield));
% %
% %Check if the new .mat has any additional fields
% if exist('cellinfo','var')    
%     matfields = fieldnames(thiscellinfo);
%     resultsfields = fieldnames(cellinfo);
%     newfields = setdiff(matfields,resultsfields);
%     if ~isempty(newfields)
%         for ff = 1:length(newfields)
%             cellinfo(1).(newfields{ff}) = []; 
%         end
%         FIELDMISMATCH=true;
%     end
% 
%     cellinfo = orderfields(cellinfo,thiscellinfo);   
% end


end

