function y = contrast(x)

if max(x(:))~=min(x(:))
    y = (x-min(x(:)))/(max(x(:))-min(x(:)));
else
    y = x./max(x(:));
end