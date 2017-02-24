function [nCoords, nDim, nVTMode, c] = plx_vt_interpret(ts, sv);
% plx_vt_interpret - interpret CinePlex video tracking data
%
% [nCoords, nDim, nVTMode, c] = plx_vt_interpret(ts, sv);
%
% INPUT:
%   ts - array of timestamps (in seconds) (see plx_event_ts.m)
%   sv - array of strobed event values (see plx_event_ts.m)
%
% OUTPUT:
%   nCoords - number of produced coordinates
%   nDim    - number of elemnts in produced coordinates
%             nDim = 3 for CENTROID, LED_1, LED_2, LED3
%             nDim = 4 for CENTROID_WITH_MOTION
%             nDim = 5 for LED_12, LED_23, LED_13
%             nDim = 7 for LED_123
%   nVTMode - VT mode:
%              0 = UNKNOWN
%              1 = CENTROID                // 1 set of coordinates, no motion
%              2 = CENTROID_WITH_MOTION    // 1 set of coordinates, with motion
%              3 = LED_1                   // 1 set of coordinates
%              4 = LED_2                  
%              5 = LED_3
%              6 = LED_12                  // 2 sets of coordinates
%              7 = LED_13
%              8 = LED_23
%              9 = LED_123                 // 3 sets of coordinates
%   c       - nCoords by nDim matrix of produced coordinates
%             c(:, 1) - timestamp
%             c(:, 2) - x1
%             c(:, 3) - y1
%             c(:, 4) - x2 or motion (if present)
%             c(:, 5) - y2 (if present)
%             c(:, 6) - x3 (if present)
%             c(:, 7) - y3 (if present)
%

if nargin < 2
    error 'Expected 2 input arguments';
end

[nCoords, nDim, nVTMode, c] = mexPlex(21, '', ts, sv);