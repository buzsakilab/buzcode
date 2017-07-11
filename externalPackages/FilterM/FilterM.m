function [Y, Z] = FilterM(b, a, X, Z, Dim, Reverse)
% FilterM - Fast C-Mex filter
% [Y, Z] = FilterM(b, a, X, Z, Dim, Reverse)
% INPUT:
%   b, a: Filter parameters as DOUBLE vectors. If the vectors have different
%      lengths, the shorter one is padded with zeros.
%   X: Signal as DOUBLE or SINGLE vector or array.
%   Z: Initial conditions as DOUBLE or SINGLE array. The size must be:
%        [(Order) - 1, SIZE(X,2), ..., SIZE(X, NDIMS(X))]
%      Optional, default: Zeros.
%   Dim: Dimension to operate on. Optional, default: 1st non-singelton.
%   Reverse: Set this to TRUE or 'reverse' to process the signal reverse order.
%      The parameters b, a and Z are not affected - see examples.
%      This is supports a faster FILTFILT operation.
%      Optional, default: FALSE.
%
% OUTPUT:
%   Y: Filtered signal with the same size and type as X.
%   Z: Final conditions as DOUBLE (!) array.
%
% NOTES:
%   - The output equals the output of FILTER exactly for DOUBLEs.
%   - For signals in SINGLE format, all intermediate values are DOUBLEs to
%     increase the accuracy. Therefore the output Z is a DOUBLE also.
%   - The parameters [a] and [b] are normalized to the 1st element of a, if
%     a(1) differs from 1.0.
%   - To accelerate my FEX submission FiltFiltM, I've implemented a filter as
%     C-Mex, which works in reverse order. To my surprise this was faster than
%     running Matlab's FILTER forward, e.g. 3.7 times for a [10000 x 1] vector, 
%     5th order Butterworth filter (Matlab 2009a, WinXP 32 bit, single core). 
%     Therefore I've expanded the Mex such that the direction can be defined
%     as input. The algorithm is a direct form II transposed structure.
%   - If the 1st dimension is processed, calling the Mex directly is faster:
%       [Y, Z] = filterX(b, a, X, Z, Reverse);
%     Z, Reverse are optional. See help section of filterX.c for instructions.
%   - In opposite to Matlab's FILTER, Z cannot be a vector if the signal is an
%     array. This is more logical.
%
% INSTALLATION:
%   Setup the compiler if not done before: mex -setup.
%   Auto-compilation: Call FilterM without inputs to start the compilation.
%   See FilterX.c for instructions for a manual compilation.
%
% EXAMPLES:
%   B = [0.000416546, 0.00124964, 0.00124964, 0.000416546];
%   A = [1, -2.68616, 2.41966, -0.730165];      % Lowpass Butterworth
%   t = linspace(0, 2*pi, 400);
%   x = sin(t) + rand(1, numel(t)) * 0.1;
%   y = filterX(B, A, x(:));     % Column vector!
%   plot(t, x, 'r', t, y, 'b');  hold('on');
%   yr = filterX(B, A, y(:), [], 'reverse');
%   plot(t, yr, 'k');
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP, 32bit
%         Compiler: LCC2.4/3.8, BCC5.5, OWC1.8, MSVC2008
% Assumed Compatibility: higher Matlab versions, Mac, Linux, 64bit
% Author: Jan Simon, Heidelberg, (C) 2011 matlab.THISYEAR(a)nMINUSsimon.de
%
% See also FILTER, FILTFILT, CONV, FiltFiltM.

% $JRev: R-g V:006 Sum:aCeGYUtdTVIU Date:20-Jul-2011 10:29:21 $
% $License: NOT_RELEASED $
% $UnitTest: uTest_FilterM $
% $File: Tools\GLMath\FilterM.m $
% History:
% 004: 19-Jul-2011 11:20, Wrapper for filterX.mex, which operates on Dim=1 only.

% Naive direct form II implementation for vectors in Matlab:
%   function [Y, z] = emulateFILTER(b, a, X, z)
%   n    = length(a);
%   z(n) = 0;
%   Y    = zeros(size(X));
%   for m = 1:length(Y)
%      Xm   = X(m);
%      Y(m) = b(1) * Xm + z(1);
%      Ym   = Y(m);
%      for i = 2:n
%         z(i - 1) = b(i) * Xm + z(i) - a(i) * Ym;
%      end
%   end
%   z = z(1:n - 1);
%
% Improvements:
% - Unroll last iteration of inner loop and omit z(n)==0.
% - Unroll inner loop completely for used filter length.

% Initialize: ==================================================================
% Global Interface: ------------------------------------------------------------
% Initial values: --------------------------------------------------------------
% Program Interface: -----------------------------------------------------------
nArg = nargin;
if nArg < 3
   if nargin == 0
      Compile;
      return;
   else
      error(['JSimon:', mfilename, ':BadNInput'], ...
         ['*** ', mfilename, ': 3 to 6 inputs required.']);
   end
elseif nArg < 4
   Z       = [];
   Dim     = [];
   Reverse = false;
elseif nArg < 5
   Dim     = [];
   Reverse = false;
elseif nArg < 6
   Reverse = false;
elseif nArg > 6
   error(['JSimon:', mfilename, ':BadNInput'], ...
      ['*** ', mfilename, ': 3 to 6 inputs required.']);
end

replyZ = (nargout >= 2);

% User Interface: --------------------------------------------------------------
% Do the work: =================================================================
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
         Z = Z(:);
      end
   elseif Dim < 3
      doReshape = false;
      doPermute = (Dim == 2);
      if doPermute
         perm = [2, 1];
         X    = transpose(X);
         Z    = transpose(Z);
      end
   else
      % Usually Matlab functions allow appending arbitrary trailing singelton
      % dimensions...
      error(['JSimon:', mfilename, ':BadDim'], ...
         ['*** ', mfilename, ': Dim out of range.']);
   end
   
else                % Signal has more than 2 dimensions: -----------------------
   if isempty(Dim)                                % 1st non-singleton dimension
      doPermute = false;
      nonSingle = find(sizeX ~= 1);
      if ~isempty(nonSingle) && nonSingle(1) > 1  % Short-circuit needed
         tmp       = [sizeX, 1];                  % RESHAPE ignores trailing 1s
         X         = reshape(X, tmp(nonSingle(1):length(tmp)));
         doReshape = true;
      else
         doReshape = false;
      end
   elseif Dim <= ndimsX                           % Dim specified:
      if Dim == 1                                 % Nothing to do on 1st dim
         doReshape = false;
         doPermute = false;
      else                                        % Specified dim to front:
         doReshape = false; % true;
         doPermute = true;
         perm      = [Dim, 1:Dim - 1, Dim + 1:ndims(X)];
         % X         = reshape(permute(X, perm), sizeX(Dim), []);
         X         = permute(X, perm);
         sizeX     = sizeX(perm);
         if ~isempty(Z)                           % Initial conditions:
            Z = permute(Z, perm);
            % Z = reshape(permute(Z, perm), size(Z, Dim), []);
         end
      end
   else
      % Usually Matlab functions allow appending arbitrary trailing singelton
      % dimensions...
      error(['JSimon:', mfilename, ':BadDim'], ...
         ['*** ', mfilename, ': Dim out of range.']);
   end
end

% Call the calculator: ---------------------------------------------------------
if replyZ
   [Y, Z] = FilterX(b, a, X, Z, Reverse);
else
   Y = FilterX(b, a, X, Z, Reverse);
end

% Reconstruct original array dimension: ----------------------------------------
if doReshape
   Y = reshape(Y, sizeX);
end
if doPermute
   Y = ipermute(Y, perm);
   if replyZ
      Z = ipermute(Z, perm);
   end
end

% return;


% ******************************************************************************
function Compile
% Compiler filterX.c

[FuncPath, FuncName] = fileparts(mfilename('fullpath'));
bakCD = cd;

% Search the C-sources:
if isempty(which(fullfile(FuncPath, 'FilterX.c')))
   error(['JSimon:', FuncName, ':Compile'], ...
      ['*** ', FuncName, ': Cannot find the C-file "FilterX.c" ???'])
end

if isunix
   % Enable the "modern" C99. Matlab defines the primordial C89 style as default
   % for the GCC compiler. If you think, that this is hilarious, send an
   % enhancement report to support@mathworks.com.
   Flags = {'-O', 'CFLAGS="\$CFLAGS', '-std=c99"', 'FilterX.c'};
else
   Flags = {'-O', 'FilterX.c'};
end

% Display the compilation command:
cmd = ['mex', sprintf(' %s', Flags{:})];
fprintf('%s\n', cmd);

try
   cd(FuncPath);
   mex(Flags{:});
   fprintf('Compilation successful.\n');
   
   % Suggest to run the unit-test:
   if ~isempty(which('uTest_FilterM'))
      if usejava('jvm')
         fprintf(['Test validity and speed: ', ...
            '<a href="matlab:uTest_FilterM">uTest_FilterM</a>\n']);
      else
         fprintf('Test validity and speed: uTest_FilterM\n');
      end
   else
      fprintf(2, '??? Cannot find unit-test uTest_FilterM\n');
   end
catch
   fprintf(2, '\n*** Compilation failed:\n%s', lasterr);
   fprintf(['\nTry to compile manually. ', ...
      'For copy&paste in the command window:\n', ...
      '  cd(''%s'')\n  %s\n'], FuncPath, cmd);
end

% Restore original directory:
cd(bakCD);

% return;
