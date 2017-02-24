function bool = CheckTS(varargin)

% bool = [c]tsd/CheckTS(X0, X1, X2, ...)
%
% checks to make sure that all timestamps are identical for all tsds included.
% works with combinations of ctsd and tsd

% ADR 1998
% version L4.0
% status: PROMOTED

adrlib;

R0 = Range(varargin{1},'ts');
for iX = 2:length(varargin)
   R1 = Range(varargin{iX},'ts');
   if (length(R0) ~= length(R1))
      bool = false;                  % if not same length, can't be equal.
      return
   end
   if (min(R0 == R1) == 0)
      bool = false;                  % if there are any non-equal elts, not equal
      return
   end
end
bool = true;                        % nothing failed, must be ok

         