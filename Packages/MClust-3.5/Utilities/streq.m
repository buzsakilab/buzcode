function bool = streq(S1, S2)

% bool = streq(S1, S2)
%
% returns true if S1 == S2 and false otherwise
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

adrlib;
bool = strcmp(S1, S2) == 1;
