function Errors = TestSelf;

% Errors = tsd/TestSelf;

adrlib;
Errors = {};

t = 1:10;
xd = reshape(1:60, 10, 2, 3);
x = tsd(t, xd);

% Data
if Data(x) ~= xd
   Errors = [Errors, {'Data(x) failed.'}];
end
if reshape(Data(x, [5 6]), 1,12) ~=  [5 6 15 16 25 26 35 36 45 46 55 56]
   Errors = [Errors, {'Data(x, []) failed.'}];
end

% CheckTS
if ~CheckTS(x,x)
   Errors = [Errors, {'CheckTS failed.'}];
end

% DT
if DT(x) ~= 1
   Errors = [Errors, {'DT failed.'}];
end

% EndTime
if EndTime(x) ~= 10
   Errors = [Errors, {'EndTime failed.'}];
end

% StartTime
if StartTime(x) ~= 1
   Errors = [Errors, {'StartTime failed.'}];
end

% Range
if Range(x,'ts') ~= t
   Errors = [Errors, {'Range failed.'}];
end

% Restrict
y = Restrict(x, [5 6]);
if reshape(Data(y),1,12) ~=  [5 6 15 16 25 26 35 36 45 46 55 56]
   Errors = [Errors {'Restrict(x, []) failed.'}];
end
y = Restrict(x, 5, 6);
if reshape(Data(y),1,12) ~=  [5 6 15 16 25 26 35 36 45 46 55 56]
   Errors = [Errors {'Restrict(x, t0, t1) failed.'}];
end

% cat
x1 = x;
x2 = ctsd(EndTime(x1)+10, 1, xd);
if Data(cat(x1,x2)) ~= cat(1, xd,xd)
   Errors = [Errors, {'cat failed.'}];
end

% mask
y = tsd(1:10,1:10);
M = Mask(y, [2 4]);
if Data(M) ~= [nan 2 3 4 nan nan nan nan nan nan]
   Errors = [Errors, {'mask failed.'}];
end

