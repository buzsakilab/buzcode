function y = contrast(x,varargin)
% 
% y = contrast(x)
% makes values of y lying in the range [0 1] from the COLUMN vector
% 
% y = contrast(x,minVal,maxVal)
% makes values of y lying in the range [minVal maxVal]. Any value outside
% those boundaries will be set to boundary value
% 
% Adrien Peyrache 2012

if size(x,1)<2
    error('Constrast needs at least 2 values')
end
if size(x,2)>1
    y = [];
    for ii=1:size(x,2)
        if ~isempty(varargin)
            y = [y contrast(x(:,ii),varargin{1},varargin{2})];
        else
            y = [y contrast(x(:,ii))];
        end
    end
else

    if isempty(varargin)
        mx = max(x(:));
        mn = min(x(:));
    elseif length(varargin)==2
        mn = varargin{1};
        mx = varargin{2};
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
end