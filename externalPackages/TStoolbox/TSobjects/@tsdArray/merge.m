function SO = merge(S1,S2,offset)

% Merges tsd's in tsdArray
% 

if length(S1)~=length(S2)
    error('must the same length')
end

SO = cell(length(S1),1);

for ii=1:length(S1)
    keyboard
    SO{ii} = ts([Range(S1.C{ii});Range(S2.C{ii})+offset]);
end
SO = tsdArray(SO);

    
  