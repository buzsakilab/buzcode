function IMR = calcImr(m,im)
%CALCIMR Calculates the Ideal Mask Ratio (IMR)
% 
%   IMR = IOSR.BSS.CALCIMR(M,IM) calculates the ideal mask ratio (IMR) from
%   the calculated mask M and ideal mask IM. Masks can be logical or
%   double, but must only contain values in the interval [0,1].

%   Copyright 2016 University of Surrey.

    % check input
    assert(nargin>1, 'iosr:calcImr:invalidInput', 'You need 2 masks to calculate IMR!')
    assert(size(m,1)>1 & size(m,2)>1, 'iosr:calcImr:invalidMask', 'm must be a two-dimensional matrix')
    assert(size(im,1)>1 & size(im,2)>1, 'iosr:calcImr:invalidMask', 'im must be a two-dimensional matrix')
    assert(all(size(m)==size(im)), 'iosr:calcImr:invalidMask', 'Masks must be the same size!')

    % validate masks
    m = check_mask(m);
    im = check_mask(im);

    % calculate IMR
    lamda = sum(sum(m.*im));
    rho = sum(sum(abs(m-im)));

    IMR = lamda/(lamda+rho);

end

function m = check_mask(m)
%CHECK_MASK validate mask

    if ~islogical(m)
        if any(m(:)<0) || any(m(:)>1);
            error('iosr:calcImr:maskValuesOutOfRange','Values outside of [0,1] were found.')
        end
    else % make numeric
        m = +m;
    end

end
