function Y = FiltFiltM(b, a, X, Dim)
% FiltFiltM - zero phase-shift digital filter
% The signal is filtered in forward and reverse direction for a zero phase
% distortion. This is equivalent to FILTFILT of Matlab's Signal Processing
% Toolbox, but the direct processing of arrays let this function run 10% to 90%
% faster. With the fast Mex function FilterX, FiltFiltM is 80% to 95% faster
% than FILTFILT and uses much less temporary memory. 
%
% Y = FiltFiltM(b, a, X, Dim)
% INPUT:
%   b, a: Filter coefficients as vectors. See FILTER.
%   X:    Signal as DOUBLE or SINGLE array.
%   Dim:  Dimension to operate on. If Dim is the empty matrix [], the filter is
%         applied to the 1st non-singelton dimension. Filtering along the 1st
%         dimension is faster than for the other dimensions.
%         Optional, default: 1.
% OUTPUT:
%   Y:    Array of same size and type as X. The values are filtered in forward
%         and backward direction for a zero phase shift.
%         The result equals the output of FILTFILT.
%
% NOTES:
% - In opposite to FILTFILT the initial and final phases are treated separately
%   instead of enlarging the signal array. This reduces the memory consumption.
% - The results for signals of type SINGLE are more accurate, because FilterX
%   uses DOUBLE values for computing the filter.
%
% EXAMPLES:
%   b = [0.0675, 0.1349, 0.0675];  % Low pass Butterworth filter
%   a = [1.0,   -1.1430, 0.4128];
%   X = rand(16, 8);
%   Y = FiltFiltM(b, a, X);        % Along the 1st dimension
%   figure, subplot(1,2,1); plot(Y);
%   Y = FiltFiltM(b, a, X, 2);     % Along the 2nd dimension
%   subplot(1,2,2), plot(Y.');
%
% INSTALLATION:
%   This function runs without installation. But for more speed the C-mex
%   FilterX.c can be compiled: Call FilterM() without inputs. See there for
%   further instructions.
%
% See also FILTER, FILTER2, FILTFILT, CONV, FilterM.
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP
% Author: Jan Simon, Heidelberg, 2009-2011 matlab.THISYEAR(a)nMINUSsimon.de

% References:
% [1] Sanjit K. Mitra, Digital Signal Processing, 2nd ed., McGraw-Hill, 2001
% [2] Fredrik Gustafsson, Determining the initial states in forward-backward
%     filtering, IEEE Transactions on Signal Processing, pp. 988--992,
%     April 1996, Volume 44, Issue 4

% $JRev: R-t V:019 Sum:flRzmua2/Vps Date:20-Jul-2011 10:29:21 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $UnitTest: uTest_FiltFiltM $
% $File: Tools\GLMath\FiltFiltM.m $
% History:
% 009: 08-May-2011 16:45, Get initial and final parts separately: 10% faster.
%      Joining the reflected initial and final parts to the signal wastes
%      memory. Now the reflected edge parts are filtered separately just to get
%      the initial conditions for filtering the original signal.
% 011: 12-Jul-2011 16:40, 4th input [Dim]. Rewritten from scratch.
%      Now this function does not contain any parts of the original FILTFILT
%      anymore. The initial conditions are determined following Gustafsson's
%      article.

% Initialize: ==================================================================
% Global Interface: ------------------------------------------------------------
% Use the much faster C-Mex replacement for FILTER:
persistent hasFilterX
if isempty(hasFilterX)
   hasFilterX = ~isempty(which('FilterX'));
   if ~hasFilterX && ~isempty(which('FilterM'));
      fprintf(['::: %s: Matlab''s FILTER is used.\n', ...
         '    Run FilterM without inputs to compile FilterX.c\n'], ...
         mfilename);
   end
end

% Initial values: --------------------------------------------------------------
% Program Interface: -----------------------------------------------------------
nArg = nargin;
if nArg == 3
   Dim = [];    % Default is specified later: 1st or 1st non-singelton dimension
elseif nArg == 4
   if numel(Dim) > 1 || ~isequal(Dim, round(abs(Dim)))
      error(['JSimon:', mfilename, ':BadDim'], ...
         ['*** ', mfilename, ': Dim must be a valid dimension.']);
   end
else
   error(['JSimon:', mfilename, ':BadNInput'], ...
      ['*** ', mfilename, ': 3 or 4 inputs required.']);
end

if isempty(b) || isempty(a) || isempty(X)
   Y = X([]);
   return;
end

% Get filter parameter vectors as row vectors and padded with zeros on demand:
nb    = numel(b);
na    = numel(a);
b     = b(:);
a     = a(:);
Order = max(nb, na);
if nb < Order
   b(Order, 1) = 0;
end
if na < Order
   a(Order, 1) = 0;
end
nEdge = 3 * (Order - 1);   % Number of edge elements

% Move dimension to operate on to dim 1: ---------------------------------------
% Specifying the dimension in the FILTER command is not faster. It seems like
% the permutation is performed inside FILTER also instead of operating on the
% original data.
sizeX  = size(X);
ndimsX = ndims(X);
if ndimsX == 2             % Vector or matrix
   if isempty(Dim)         % Default: First non-singelton dimension
      doPermute = false;
      doReshape = (sizeX(1) == 1);
      if doReshape
         X = X(:);
      end
   elseif Dim < 3
      doReshape = false;
      doPermute = (Dim == 2);
      if doPermute
         perm = [2, 1];
         X    = transpose(X);
      end
   else
      % Usually Matlab functions allow appending arbitrary trailing singelton
      % dimensions, but then the (xLen<=nEdge) check fails anyway:
      error(['JSimon:', mfilename, ':BadDim'], ...
         ['*** ', mfilename, ': Dim out of range.']);
   end
   
else                % Signal has more than 2 dimensions: -----------------------
   doReshape = true;
   if isempty(Dim)         % First non-singleton dimension:
      doPermute = false;
      nonSingle = find(sizeX ~= 1);
      if ~isempty(nonSingle) && nonSingle(1) > 1  % Short-circuit needed
         X = reshape(X, sizeX(nonSingle(1)), []);
      end
   elseif Dim <= ndimsX    % Dim specified:
      if Dim == 1
         doPermute = false;
      else
         doPermute = true;
         perm      = [Dim, 1:Dim - 1, Dim + 1:ndims(X)];
         X         = reshape(permute(X, perm), sizeX(Dim), []);
         sizeX     = sizeX(perm);
      end
   else
      error(['JSimon:', mfilename, ':BadDim'], ...
         ['*** ', mfilename, ': Dim out of range.']);
   end
end

xLen = size(X, 1);          % Not sizeX - would fail for row vectors!
if xLen <= nEdge            % Input data too short
   error(['JSimon:', mfilename, ':BadDataLen'], ...
      ['*** ', mfilename, ...
      ': Signal must be longer than 3 times the filter order.']);
end

% User Interface: --------------------------------------------------------------
% Do the work: =================================================================
% Create initial conditions to treat offsets at beginning and end:
K       = eye(Order - 1);
K(:, 1) = a(2:Order);
K(1)    = K(1) + 1;
K(Order:Order:numel(K)) = -1.0;
IC      = K \ (b(2:Order) - a(2:Order) * b(1));

% Use a reflection to extrapolate signal at beginning and end to reduce edge
% effects (BSXFUN would be some micro-seconds faster, but it is not available in
% Matlab 6.5):
x1_2 = 2 * X(1,    :);
xf_2 = 2 * X(xLen, :);
Xi   = x1_2(ones(1, nEdge), :) - X((nEdge + 1):-1:2, :);
Xf   = xf_2(ones(1, nEdge), :) - X((xLen - 1):-1:(xLen - nEdge), :);

if hasFilterX  % Use the faster C-mex filter function: -------------------------
   % Filter initial reflected signal:
   [dum, Zi] = FilterX(b, a, Xi, IC * Xi(1, :));   %#ok<ASGLU>
   
   % Use the final conditions of the initial part for the actual signal:
   [Ys, Zs]  = FilterX(b, a, X,  Zi);              % "s"teady state
   Yf        = FilterX(b, a, Xf, Zs);              % "f"inal conditions
   
   % Filter signal again in reverse order:
   [dum, Zf] = FilterX(b, a, Yf, IC * Yf(nEdge, :), 'reverse');  %#ok<ASGLU>
   Y         = FilterX(b, a, Ys, Zf, 'reverse');
   
else  % Use the slower built-in FILTER function: -------------------------------
   [dum, Zi] = filter(b, a, Xi, IC * Xi(1, :));     %#ok<ASGLU>
   [Ys, Zs]  = filter(b, a, X,  Zi);                % "s"teady state
   Yf        = filter(b, a, Xf, Zs);                % "f"inal conditions
   Yf        = Yf(nEdge:-1:1, :);
   [dum, Zf] = filter(b, a, Yf, IC * Yf(1, :));     %#ok<ASGLU>
   Y         = filter(b, a, Ys(xLen:-1:1, :), Zf);  % Filter reverted signal
   Y         = Y(xLen:-1:1, :);                     % Re-revert
end

% Reconstruct original array dimension: ----------------------------------------
if doReshape
   Y = reshape(Y, sizeX);
end
if doPermute
   Y = ipermute(Y, perm);
end

% return;
