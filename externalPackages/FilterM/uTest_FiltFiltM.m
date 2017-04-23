function uTest_FiltFiltM(doSpeed)
% Unit test: FiltFiltM
% This is a routine for automatic testing. It is not needed for processing and
% can be deleted or moved to a folder, where it does not bother.
%
% uTest_FiltFiltM(doSpeed)
% INPUT:
%   doSpeed: If ANY(doSpeed) is TRUE, the speed is measured. Optional.
% OUTPUT:
%   On failure the test stops with an error.
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP
% Author: Jan Simon, Heidelberg, (C) 2009-2011 matlab.THISYEAR(a)nMINUSsimon.de

% $JRev: R-k V:010 Sum:VKPp090JCPL7 Date:19-Jul-2011 19:40:07 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $File: Tools\UnitTests_\uTest_FiltFiltM.m $

% Init: ========================================================================
ErrHead = ['JSimon:', mfilename, ':TestFailed'];

if nargin == 0 || any(doSpeed)
   TestTime = 2.0;  % sec
else
   TestTime = 0.5;  % sec
end

% ==============================================================================
disp(['==== Test FiltFiltM  ', datestr(now, 0)]);
disp(['  Version: ', which('FiltFiltM')]);
pause(0.01);  % Flush events

fprintf('\nKnown answer tests:\n');
B = [];
A = [];
Y = FiltFiltM(B, A, []);
if ~isempty(Y)
   error(ErrHead, ['*** ', mfilename, ': Bad reply for empty signal']);
end
disp('  ok: all input empty');

B = [];
A = [];
X = rand(1000, 1);
Y = FiltFiltM(B, A, X);
if ~isempty(Y)
   error(ErrHead, ['*** ', mfilename, ': Bad reply for empty parameters']);
end
disp('  ok: empty parameters');

% [B, A] = ButterParam(2, 0.2, 'low')
B = [0.06745527, 0.1349105, 0.06745527];
A = [1, -1.142981, 0.4128016];
Y = FiltFiltM(B, A, []);
if ~isempty(Y)
   error(ErrHead, ...
      ['*** ', mfilename, ': Bad reply for empty signal with parameters']);
end
disp('  ok: empty signal');

Y = FiltFiltM(B, A, X);
Z = filtfilt(B, A, X);
if ~isequal(size(Y), size(Z)) || any(abs(Y - Z) > 10 * eps)
   error(ErrHead, ...
      ['*** ', mfilename, ': Bad result for [1000 x 1]']);
end
disp('  ok: [N x 1] vector');

Xt = reshape(X, 1, []);
Y  = FiltFiltM(B, A, Xt);
if ~isequal(size(Y), size(Xt)) || any(abs(Y(:) - Z(:)) > 10 * eps)
   error(ErrHead, ...
      ['*** ', mfilename, ': Bad result for [1 x 1000]']);
end
disp('  ok: [1 x N] vector');

X = [X, X];
Y = FiltFiltM(B, A, X);
Z = filtfilt(B, A, X);
if ~isequal(size(Y), size(Z)) || any(abs(Y(:) - Z(:)) > 10 * eps)
   error(ErrHead, ...
      ['*** ', mfilename, ': Bad result for [1000 x 2]']);
end
disp('  ok: 2D array, 1st dim');

X = transpose(X);
Y = FiltFiltM(B, A, X, 2);
Z = transpose(Z);
if ~isequal(size(Y), size(Z)) || any(abs(Y(:) - Z(:)) > 10 * eps)
   error(ErrHead, ...
      ['*** ', mfilename, ': Bad result for [2 x 1000], 2nd dim']);
end
disp('  ok: 2D array, 2nd dim');

X = rand(20, 30, 40);
Y = FiltFiltM(B, A, X);
Z = filtfilt(B, A, X);
if ~isequal(size(Y), size(Z)) || any(abs(Y(:) - Z(:)) > 10 * eps)
   error(ErrHead, ...
      ['*** ', mfilename, ': Bad result for [20 x 30 x 40]']);
end

Y = FiltFiltM(B, A, X, 2);
Z = ipermute(filtfilt(B, A, permute(X, [2,1,3])), [2,1,3]);
if ~isequal(size(Y), size(Z)) || any(abs(Y(:) - Z(:)) > 10 * eps)
   error(ErrHead, ...
      ['*** ', mfilename, ': Bad result for [20 x 30 x 40], 2nd dim']);
end

Y = FiltFiltM(B, A, X, 3);
Z = ipermute(filtfilt(B, A, permute(X, [3,1,2])), [3,1,2]);
if ~isequal(size(Y), size(Z)) || any(abs(Y(:) - Z(:)) > 10 * eps)
   error(ErrHead, ...
      ['*** ', mfilename, ': Bad result for [20 x 30 x 40], 3rd dim']);
end
disp('  ok: 3D array, 1st, 2nd, 3rd dim');

X = rand(1000, 1);
Y = FiltFiltM(B, A, X);
Z = FiltFiltM(B, A, reshape(X, 1, 1, []), []);
if ~isequal(size(Z), [1, 1, numel(X)]) || any(abs(Y(:) - Z(:)) > 10 * eps)
   error(ErrHead, ['*** ', mfilename, ...
      ': Bad result for [1 x 1 x 1000], 1st non-singelton']);
end
disp('  ok: [1 x 1 x N] array, 1st non-singelton');

% FILTER does not support SINGLES in Matlab 6.5:
try
   y = filter(B, A, single(ones(1, 100)));   %#ok<NASGU>
   FILTERforSINGLE = true;
catch
   FILTERforSINGLE = false;
   lasterr('');
end

if FILTERforSINGLE
   X  = rand(1000, 2);
   Y  = FiltFiltM(B, A, X);
   Xs = single(X);
   Ys = FiltFiltM(B, A, Xs);
   Zs = filtfilt(B, A, Xs);
   if ~isequal(size(Ys), size(Zs))
      error(ErrHead, ['*** ', mfilename, ': Bad result for SINGLE']);
   end
   epsSingle    = 1.192092895507813e-007;
   errFiltFiltM = max(abs(double(Ys(:)) - Y(:))) / epsSingle;
   errFILTFILT  = max(abs(double(Zs(:)) - Y(:))) / epsSingle;
   if any(errFiltFiltM > 10)
      error(ErrHead, ['*** ', mfilename, ': Bad result for SINGLE']);
   end
   disp('  ok: SINGLE array');
   fprintf(['      Difference to DOUBLE divided by EPS(SINGLE(1)):\n', ...
      '      FilterM: %g, FILTER: %g\n'], errFiltFiltM, errFILTFILT);
   
else
   v = sscanf(version, '%d.', 2);
   fprintf('  ok: FILTER does not accept SINGLE in Matlab %d.%d\n', v);
end

% ------------------------------------------------------------------------------
tooLazy = false;
fprintf('\nProvoke problems:\n');
try
   Y = FiltFiltM(B, A, rand(4));  %#ok<NASGU>
   tooLazy = true;
catch
   disp('  ok: Too short signal in dim 1 caught.');
end
if tooLazy
   error(ErrHead, ['*** ', mfilename, ': Too short signal accepted.']);
end

try
   Y = FiltFiltM(B, A, rand(4), 2);  %#ok<NASGU>
   tooLazy = true;
catch
   disp('  ok: Too short signal in dim 2 caught.');
end
if tooLazy
   error(ErrHead, ['*** ', mfilename, ': Too short signal accepted.']);
end

try
   Y = FiltFiltM(B, A, rand(1, 100), 3);  %#ok<NASGU>
   tooLazy = true;
catch
   disp('  ok: Dim out of range caught.');
end
if tooLazy
   error(ErrHead, ['*** ', mfilename, ': Dim out of range accepted.']);
end

try
   Y = FiltFiltM(B, A, rand(1, 100), -3.1);  %#ok<NASGU>
   tooLazy = true;
catch
   disp('  ok: Not integer Dim caught.');
end
if tooLazy
   error(ErrHead, ['*** ', mfilename, ': Not integer Dim accepted.']);
end

% Time test: -------------------------------------------------------------------
% Results on a Pentium-M 1500 MHz, 512MB RAM, WinXP, Matlab 2009a:
% Number of loops to get about 2 seconds run time
% 1000000 x 1:  6 loops,  FILTFILT: 1.930 sec,  FiltFiltM: 0.385 sec  ==> 20.0%
% 100000 x 1:  70 loops,  FILTFILT: 2.077 sec,  FiltFiltM: 0.488 sec  ==> 23.5%
% 10000 x 10:  82 loops,  FILTFILT: 1.996 sec,  FiltFiltM: 0.600 sec  ==> 30.0%
% 1000 x 100:  50 loops,  FILTFILT: 2.761 sec,  FiltFiltM: 0.518 sec  ==> 18.7%
% 100 x 100:  100 loops,  FILTFILT: 1.935 sec,  FiltFiltM: 0.079 sec  ==>  4.1%
% 100 x 1000:  12 loops,  FILTFILT: 2.303 sec,  FiltFiltM: 0.130 sec  ==>  5.6%

% A 4th order Butterworth filter: [B, A] = Butter(4, 0.25, 'low')
B = [0.0102094807912032, 0.0408379231648126, 0.0612568847472189, ...
   0.0408379231648126, 0.0102094807912032];
A = [1.0, -1.96842778693852, 1.73586070920889, -0.724470829507363, ...
   0.120389599896245];
fprintf('\n== Speed test with Butterworth filter of order %d:\n', ...
   length(B) - 1);

fprintf('Number of loops determined for about %g seconds run time\n', TestTime);

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
      origy = filtfilt(B, A, x);  %#ok<NASGU>
      clear('origy');
   end
   nLoop = max(2, ceil(nLoop * TestTime));
   fprintf('%5d loops\n', nLoop);
   
   % Run test for vectorized function:
   tic;
   for i = 1:nLoop
      origy = filtfilt(B, A, x);  %#ok<NASGU>
      clear('origy');
   end
   elapsed0 = toc;
   fprintf('    FILTFILT:  %8.3f sec\n', elapsed0);
   
   % Run test for vectorized function:
   tic;
   for i = 1:nLoop
      newy = FiltFiltM(B, A, x);  %#ok<NASGU>
      clear('newy');
   end
   elapsed1 = toc;
   fprintf('    FiltFiltM: %8.3f sec', elapsed1);
   fprintf('  ==> %.1f%%\n', 100 * elapsed1 / elapsed0);
   
   % Compare the results - usually the results are identical, but small round
   % off errors would be ok:
   origy = filtfilt(B, A, x);
   newy  = FiltFiltM(B, A, x);
   if isequal(origy, newy)
      fprintf('    Results are identical\n\n');
   elseif all(abs(origy(:) - newy(:)) < 10 * eps)
      fprintf('    Results are equal within 10*eps\n\n');
   else
      error(ErrHead, 'Different results - DO NOT USE FiltFiltM!');
   end
end

% Bye:
disp('== FiltFiltM passed the tests.');

return;
