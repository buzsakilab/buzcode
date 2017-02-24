function [yu xu] = upsample(y,n,x)

% Upsample signa

n=round(n);
l = length(y);

yu = zeros(l*n,1);

if size(y,1) == 1
	y = y';
end

yu(1:n:end) = y;
xu = [x(1):(x(end)-x(1))/(length(yu)-1):x(end)];

d = median(diff(xu));

sc = sinc(pi/n*xu/d);
yu = convn(yu,sc','same');

keyboard