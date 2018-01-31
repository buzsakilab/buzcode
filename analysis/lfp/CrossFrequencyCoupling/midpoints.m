function mids = midpoints(x)
% mids = midpoints(x)
%
% x must be a vector.

if size(x,1) == 1
    mids = zeros(1,length(x)-1);
else
    mids = zeros(length(x)-1,1);
end
x = x(:);

mids(:) = mean([x(1:(end-1)) x(2:end)], 2);
