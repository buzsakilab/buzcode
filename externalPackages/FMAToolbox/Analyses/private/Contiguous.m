%Contiguous - Determine a patch of contiguous points.
%
%  Given a matrix of zeros and ones, where ones form 'patches' of contiguous
%  items, find all items that belong to the same patch as item (i,j).
%
%  USAGE
%
%    y = Contiguous(x,i,j)
%
%    x              input matrix
%    i,j            starting coordinates (row first)
%
%  OUTPUT
%
%    Y              output matrix, where for all k and l
%                    * y(k,l) = 0 if x(k,l) = 0
%                    * y(k,l) = 1 if x(k,l) = 1 but (k,l) is not in the same
%                      patch as (i,j)
%                    * y(k,l) = 2 if x(k,l) = 1 and (k,l) is in the same patch
%                      as (i,j)

% Copyright (C) 2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
