function R = Restrict_int(tsa, t0, t1)

% tsd/Restrict
% 	R = Restrict(tsa, t0, t1)
% 	Returns a new tsa R so that tsa.data is between timestamps t0 and t1
%

% ADR 
% version L4.1
% status: PROMOTED

% v4.1 29 oct 1998 now can handle nargin=2
% v4.2 23 jan 1999 now can handle t0 and t1 as arrays

switch nargin
case 2                             % R = Restrict(tsd, t)
   error('this function is for intervals!!!');
   
case 3                             % R = Restrict(tsd, t0, t1)
   if length(t0) ~= length(t1)
      error('t0 and t1 must be same length')
   end
   ix = Restrict_intervals(tsa.t, t0, t1);
   
otherwise
   error('Unknown number of input arguments.');
   
end % switch

R = tsd(tsa.t(ix), SelectAlongFirstDimension(tsa.data, ix));
