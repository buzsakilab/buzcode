function uTest_FilterM(doSpeed)
% Unit test: FilterM
% This is a routine for automatic testing. It is not needed for processing and
% can be deleted or moved to a folder, where it does not bother.
%
% uTest_FilterM(doSpeed)
% INPUT:
%   doSpeed: If ANY(doSpeed) is TRUE, the speed is measured. Optional.
% OUTPUT:
%   On failure the test stops with an error.
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP
% Author: Jan Simon, Heidelberg, (C) 2009-2011 matlab.THISYEAR(a)nMINUSsimon.de

% *** NOTE:
% The results of MATLAB's FILTER and FilterX.mex can differ due to rounding.
% But at least with MSVC 2008 the results are identical.
% ***

% $JRev: R-m V:012 Sum:Uwvd0UKv8byD Date:20-Jul-2011 10:32:34 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $File: Tools\UnitTests_\uTest_FilterM.m $

% Init: ========================================================================
ErrHead = ['JSimon:', mfilename, ':TestFailed'];

if nargin == 0 || any(doSpeed)
   TestTime = 2.0;  % sec
else
   TestTime = 0.5;  % sec
end

% Run tests for a specific filter, implemented orders: 2,4,5,8
[B, A] = localButter(4);

% ==============================================================================
disp(['==== Test FilterM  ', datestr(now, 0)]);
disp(['  Version: ', which('FilterM')]);
whichX = which('FilterX');
if isempty(whichX)
   error(ErrHead, 'Compile the Mex at first: Call FilterM without inputs!');
end
disp(['           ', whichX]);

pause(0.01);  % Flush events

fprintf('\nKnown answer tests:\n');
Y = FilterM([], [], []);
if ~isempty(Y)
   error(ErrHead, ['*** ', mfilename, ': Bad reply for empty signal']);
end
disp('  ok: all input empty');

X = rand(1000, 1);
Y = FilterM([], [], X);
if ~isempty(Y)
   error(ErrHead, ['*** ', mfilename, ': Bad reply for empty parameters']);
end
disp('  ok: empty parameters');

Y = FilterM(B, A, []);
if ~isempty(Y)
   error(ErrHead, ...
      ['*** ', mfilename, ': Bad reply for empty signal with parameters']);
end
disp('  ok: empty signal');

Y = FilterM(B, A, X);
Z = filter(B, A, X);
if ~isequal(size(Y), size(Z)) || ~isequalTol(Y, Z)
   error(ErrHead, ...
      ['*** ', mfilename, ': Bad result for [1000 x 1]']);
end
disp('  ok: [N x 1] vector');

Xt = reshape(X, 1, []);
Y  = FilterM(B, A, Xt);
if ~isequal(size(Y), size(Xt)) || ~isequalTol(Y, transpose(Z))
   error(ErrHead, ...
      ['*** ', mfilename, ': Bad result for [1 x 1000]']);
end
disp('  ok: [1 x N] vector');

X = [X, X];
Y = FilterM(B, A, X);
Z = filter(B, A, X);
if ~isequal(size(Y), size(Z)) || ~isequalTol(Y, Z)
   error(ErrHead, ...
      ['*** ', mfilename, ': Bad result for [1000 x 2]']);
end
disp('  ok: 2D array, 1st dim');

Y2 = transpose(FilterM(B, A, transpose(X), [], 2));
if ~isequal(Y, Y2)
   error(ErrHead, ...
      ['*** ', mfilename, ': Bad result for [2 x 1000], 2nd dim']);
end
disp('  ok: 2D array, 2nd dim');

X = rand(20, 30, 40);
Y = FilterM(B, A, X);
Z = filter(B, A, X);
if ~isequal(size(Y), size(Z)) || ~isequalTol(Y, Z)
   error(ErrHead, ...
      ['*** ', mfilename, ': Bad result for [20 x 30 x 40]']);
end

Y = FilterM(B, A, X, [], 2);
Z = ipermute(filter(B, A, permute(X, [2,1,3])), [2,1,3]);
if ~isequal(size(Y), size(Z)) || ~isequalTol(Y, Z)
   error(ErrHead, ...
      ['*** ', mfilename, ': Bad result for [20 x 30 x 40], 2nd dim']);
end

Y = FilterM(B, A, X, [], 3);
Z = ipermute(filter(B, A, permute(X, [3,1,2])), [3,1,2]);
if ~isequal(size(Y), size(Z)) || ~isequalTol(Y, Z)
   error(ErrHead, ...
      ['*** ', mfilename, ': Bad result for [20 x 30 x 40], 3rd dim']);
end
disp('  ok: 3D array, 1st, 2nd, 3rd dim');

X = rand(1000, 1);
Y = FilterM(B, A, X);
Z = FilterM(B, A, reshape(X, 1, 1, []), []);
if ~isequal(size(Z), [1, 1, numel(X)]) || ~isequalTol(Y(:), Z(:))
   error(ErrHead, ['*** ', mfilename, ...
      ': Bad result for [1 x 1 x 1000], 1st non-singelton']);
end
disp('  ok: [1 x 1 x N] array, 1st non-singelton');

X  = rand(1000, 4);
Y  = FilterM(B, A, X, [], [], 'forward');
Y2 = flipud(FilterM(B, A, flipud(X), [], [], 'reverse'));
if ~isequal(Y2, Y)
   error(ErrHead, ['*** ', mfilename, ': Bad result for reverse mode']);
end
disp('  ok: reverse');

try
   y = filter(B, A, single(ones(1, 100)));   %#ok<NASGU>
   FILTERforSINGLE = true;
catch
   FILTERforSINGLE = false;
   lasterr('');
end

if FILTERforSINGLE
   Xs = single(rand(1000, 2));
   Xd = double(Xs);             % DOUBLE precision with SINGLE accuracy
   Yd = FilterM(B, A, Xd);
   Ys = FilterM(B, A, Xs);
   Zs = filter(B, A, Xs);
   if ~isequal(size(Ys), size(Zs))
      error(ErrHead, ['*** ', mfilename, ': Bad result for SINGLE']);
   end
   % eps(single(1)) is not available in Matlab 6.5:
   epsSingle  = 1.192092895507813e-007;
   errFilterM = max(abs(double(Ys(:)) - Yd(:))) / epsSingle;
   errFILTER  = max(abs(double(Zs(:)) - Yd(:))) / epsSingle;
   if any(errFilterM > 10)
      error(ErrHead, ['*** ', mfilename, ': Bad result for SINGLE']);
   end
   disp('  ok: SINGLE array');
   
   if errFILTER > 100
      Wow = ' (!!!)';
   else
      Wow = '';
   end
   fprintf(['      Difference to DOUBLE divided by EPS(SINGLE(1)):\n', ...
      '      FilterM: %g, FILTER: %g%s\n'], errFilterM, errFILTER, Wow);
   
   Ys2 = flipud(FilterM(B, A, flipud(Xs), [], [], 'reverse'));
   if ~isequal(Ys2, Ys)
      error(ErrHead, ['*** ', mfilename, ': Bad result for SINGLE, reverse']);
   end
   disp('  ok: SINGLE array, reverse');

else
   v = sscanf(version, '%d.', 2);
   fprintf('  ok: FILTER does not accept SINGLE in Matlab %d.%d\n', v);
end

% Initial conditions: ----------------------------------------------------------
identical = true;

x      = rand(20, 30, 40, 45);
order  = length(B);
zi     = rand(order - 1, 30, 40, 45);
[y, z] = FilterM(B, A, x, zi, 1);

for d = 1:4
   % Permute 1st to dimension d'th dimension, filter along d'th dimension:
   p    = 1:4;
   p(d) = 1;
   p(1) = d;
   [yp, zp] = FilterM(B, A, permute(x, p), permute(zi, p), d);
   y2 = ipermute(yp, p);
   z2 = ipermute(zp, p);
   if not(isequal(y, y2) && isequal(z, z2))
      identical = false;
      if not(isequalTol(y, y2) && isequalTol(z, z2))
         error(ErrHead, ...
            ['*** ', mfilename, ': Working on permutated array failed.']);
      end
   end
end

if identical
   disp('  ok: permutated array');
else
   disp('  ok: permutated array, round-off');  % Not expected!
end

% ------------------------------------------------------------------------------
tooLazy = false;
fprintf('\nProvoke problems:\n');
try
   Y = FilterM(B, A, rand(4), 2);  %#ok<NASGU>
   tooLazy = true;
catch
   disp('  ok: Too short signal in dim 2 caught.');
end
if tooLazy
   error(ErrHead, ['*** ', mfilename, ': Too short signal accepted.']);
end

try
   Y = FilterM(B, A, rand(1, 100), 3);  %#ok<NASGU>
   tooLazy = true;
catch
   disp('  ok: Dim out of range caught.');
end
if tooLazy
   error(ErrHead, ['*** ', mfilename, ': Dim out of range accepted.']);
end

try
   Y = FilterM(B, A, rand(1, 100), -3.1);  %#ok<NASGU>
   tooLazy = true;
catch
   disp('  ok: Not integer Dim caught.');
end
if tooLazy
   error(ErrHead, ['*** ', mfilename, ': Not integer Dim accepted.']);
end

% Time test: -------------------------------------------------------------------
% Results on a Pentium-M 1500 MHz, 512MB RAM, WinXP, Matlab 2009a:
% Loops determined for about 2 seconds run time
%   1000000 x 1:   26 loops
%     filter:     2.030 sec
%     FilterM:    0.752 sec  ==> 37.0%
%   100000 x 1:   250 loops
%     filter:     1.988 sec
%     FilterM:    0.793 sec  ==> 39.9%
%   10000 x 10:   250 loops
%     filter:     1.986 sec
%     FilterM:    0.804 sec  ==> 40.5%
%   1000 x 100:   250 loops
%     filter:     1.999 sec
%     FilterM:    0.769 sec  ==> 38.5%
%   100 x 100:   2912 loops
%     filter:     1.968 sec
%     FilterM:    0.640 sec  ==> 32.5%
%   100 x 1000:   248 loops
%     filter:     2.002 sec
%     FilterM:    0.764 sec  ==> 38.2%

% A 4th order Butterworth filter, implemented orders: 2,4,5,8
[B, A] = localButter(4);

fprintf('\n== Speed test with Butterworth filter of order %d:\n', ...
   length(B) - 1);
fprintf('Loops determined for about %g seconds run time\n', TestTime);

TestSize = {[1e6, 1], [1e5, 1], [1e4, 10], [1e3, 100], [100, 100], [100, 1e3]};
for iTest = 1:length(TestSize)
   % Create test data:
   aSize = TestSize{iTest};
   x     = rand(aSize(1), aSize(2));
   msg   = sprintf('  %d x %d:', aSize);
   fprintf('%-14s', msg);
   drawnow;
   
   % Run test of original function:
   iTime = cputime;
   nLoop = 0;
   while cputime - iTime < 1
      nLoop = nLoop + 1;
      origy = filter(B, A, x);  %#ok<NASGU>
      clear('origy');
   end
   nLoop = max(2, ceil(nLoop * TestTime));
   fprintf('%5d loops\n', nLoop);
   
   % Run test for vectorized function:
   tic;
   for i = 1:nLoop
      origy = filter(B, A, x);  %#ok<NASGU>
      clear('origy');
   end
   elapsed0 = toc;
   fprintf('    filter:  %8.3f sec\n', elapsed0);
   
   % Run test for vectorized function:
   tic;
   for i = 1:nLoop
      newy = FilterM(B, A, x);  %#ok<NASGU>
      clear('newy');
   end
   elapsed1 = toc;
   fprintf('    FilterM: %8.3f sec', elapsed1);
   fprintf('  ==> %.1f%%\n', 100 * elapsed1 / elapsed0);
   
   % Compare the results - usually the results are identical, but small round
   % off errors would be ok:
   origy = filter(B, A, x);
   newy  = FilterM(B, A, x);
   if isequal(origy, newy)
      fprintf('    Results are identical\n\n');
   elseif isequalTol(x, y)
      fprintf('    Results are equal within 10*eps\n\n');
   else
      error(ErrHead, 'Different results - DO NOT USE FilterM!');
   end
end

% Bye:
disp('== FilterM passed the tests.');

return;

% ******************************************************************************
function [B, A] = localButter(Order)
% Some arbitrary filter parameters - keep independence from SPT!

switch Order
   case 2  % [B, A] = butter(2, 0.2, 'low')
      B = [0.06745527, 0.1349105, 0.06745527];
      A = [1.0, -1.142981, 0.4128016];
   case 4  % [B, A] = butter(4, 0.25, 'low')
      B = [0.0102094807912032, 0.0408379231648126, 0.0612568847472189, ...
         0.0408379231648126, 0.0102094807912032];
      A = [1.0, -1.96842778693852, 1.73586070920889, -0.724470829507363, ...
         0.120389599896245];
   case 5  % [B, A] = butter(5, [0.8, 0.93], 'bandpass')
      B = [0.000195572470979841, 0, -0.0009778623548992052, 0, ...
         0.00195572470979841, 0, -0.00195572470979841, 0, ...
         0.0009778623548992052, 0, -0.000195572470979841];
      A = [1.0, 8.079688736325432, 29.90860709527351, 66.78177081285723, ...
         99.59422167506651, 103.6493308964164, 76.23423011159095, ...
         39.13225635802741, 13.41992798700971, 2.777561086068105, ...
         0.2636472042716279];
   case 8  % [B, A] = Butter(8, 0.93, 'high'):
      B = [1.254203230161011e-008, -1.003362584128809e-007, ...
         3.51176904445083e-007,    -7.02353808890166e-007, ...
         8.779422611127075e-007,   -7.02353808890166e-007, ...
         3.51176904445083e-007,    -1.003362584128809e-007, ...
         1.254203230161011e-008];
      A = [1.0, 6.873011023619378,  20.73659649657562, 35.86467634824609, ...
         38.88439453756806, 27.05783359216708, 11.79932960154672, ...
         2.947769270203211, 0.3229728093056353];
   otherwise
      error(['JSimon:', mfilename, ':BadSwitch'], 'Unknown test data.');
end

% return;

% ******************************************************************************
function Equal = isequalTol(x, y)
% Check equality with a small absolute distance

if isequal(x, y)  % Try the exact and fast comparison at first
   Equal = true;
   return;
end

Tol   = 10 * eps;
Equal = false;
if isequal(size(x), size(y))
   xMy = double(x) - double(y);   % No operations on SINGLEs in Matlab 6!
   
   % Same as "if all(abs(xMy(:)) <= Tol)", but faster:
   if all(or((abs(xMy) <= Tol), (x == y)))   % is FALSE for NaNs
       Equal = true;
   end
end

% return;
