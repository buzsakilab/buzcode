function [interped] = makelength(varargin)
% figure it out.
signal = varargin{1};
len = varargin{2};
if nargin > 2
    method = varargin{3};
else
    method = 'linear';
end
    interped = interp1(1:length(signal),signal,1:length(signal)/len:length(signal),method);
c=1;
while length(interped) < len
    
    interped = interp1(1:length(signal),signal,1:length(signal)/(len+c):length(signal),method);
    
    c = 1+c;
end
end