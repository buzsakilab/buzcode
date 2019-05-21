function [new_item,new_item2, order] = sort_cells(varargin)
% Sorts two matrices by ordering the maximums of the first and using that
% order to rearrange both
% #1 = max, #2 = min


if nargin == 3
    item=varargin{1};
    item2=varargin{2};
    num=varargin{3};
elseif nargin == 1
    item=varargin{1};
    item2=varargin{1};
    num=1;
end
d = zeros(size(item,1),1);
new_item = zeros(size(item,1),size(item,2));
new_item2 = zeros(size(item2,1),size(item2,2));

if num == 1
for i = 1:size(item,1)
[blah d(i)] = max(item(i,:));    
end

[dd ddd] = sort(d);

for i = 1:size(item,1)
   new_item(i,:) = item(ddd(i),:); 
   new_item2(i,:) = item2(ddd(i),:);
end

order = ddd;
end

if num == 2
        for i = 1:size(item,1)
[blah d(i)] = min(item(i,:));    
        end

[dd ddd] = sort(d);

for i = 1:size(item,1)
   new_item(i,:) = item(ddd(i),:); 
   new_item2(i,:) = item2(ddd(i),:);
end

order = ddd;

end



