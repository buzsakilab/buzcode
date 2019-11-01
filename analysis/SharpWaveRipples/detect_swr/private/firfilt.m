% FIRFILT       FIR filtering with zero phase distortion.
%
%               matrix columns are filtered;
%               output is a column vector / matrix.

% 19-Jul-02 ES
% 27-jan-03 zero phase lag

function Y = firfilt(x,W)

if all(size(W)>=2), error('window must be a vector'), end
if numel(x)==max(size(x)), x=x(:); end

C = length(W);
if C > size( x, 1 )
    Y = NaN * ones( size( x ) );
    return
end
D = ceil(C/2) - 1;
Y = filter(W,1,[flipud(x(1:C,:)); x; flipud(x(end-C+1:end,:))]);
clear x
Y = Y(1+C+D:end-C+D,:);

return

% EOF
