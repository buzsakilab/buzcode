function y = contrast2(x,varargin)
% 
% y = contrast(x)
% makes values of matrix y lying in the range [0 1] from the input matrix x
% 
% y = contrast(x,maxVal,minVal)
% makes values of y lying in the range [minVal maxVal]. Any value outside
% those boundaries will be set to boundary value
% 
% Adrien Peyrache 2012

if size(x,1)<2
    error('Constrast needs at least 2 values')
end
if isempty(varargin)
    mx = max(x(:));
    mn = min(x(:));
elseif length(varargin)==2
    mx = varargin{1};
    mn = varargin{2};
else
    error('problem with your agument list')
end

x(x>mx)=mx;
x(x<mn)=mn;

if mx~=mn
    y = (x-mn)/(mx-mn);
else
    y = x./mx;
end
