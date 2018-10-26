function [data] = bz_shuffleCircular(data)
% circularly shifts each row

for i=1:size(data,1)
    data(i,:) = circshift(data(i,:),randi(size(data,2)));
end