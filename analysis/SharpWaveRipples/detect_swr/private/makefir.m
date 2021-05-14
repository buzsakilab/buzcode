% MAKEFIR           design an n-th order fir filter using firls
%
% H = MAKEFIR( FREQ, FS, N, BAND, TOL, X, FFUNC )
%
% FREQ - [ highpass lowpass ] frequencies (Hz)
% FS - sampling frequency (Hz); {1000}
%
% optional:
% N - order (number of samples); determined automatically if not specified
% BAND - 'low', 'high', 'bandpass', 'bandstop'; determined automatically by FREQ if not specified
% TOL - tolerance for sidebands; {0.15}
%
% to also filter some data, specify:
% X - data to be filtered
% FFUNC - filter function; {FIRFILT} and FILTFILT supported
% 
% conventions for setting FREQ/BAND:
% 100 Hz lowpass filter: [ ## 100 ] (## is ignored) 
% 10 Hz highpass: [ 10 ## ] (## is ignored)
% 10-100 Hz bandpass: [ 10 100 ] (specify BAND as 'bandstop' if desired)
%
% examples:
% h = makefir( [ 4 6 ], 1250, [], 'bandstop' ); fvtool( h ), figure( 2 ), plot( h )
% will generate a 5 Hz notch filter
%
% see also: FIRFILT

% 01-jun-12 ES

% NOTE (1) that this generates a filter with an even number of samples (see
% firls for a getaround)
% (2) this generates non-unity firs

function [ h, xf ] = makefir( freq, Fs, n, band, tol, x, ffunc )

% constants (for filter order)
FACTOR = 3;
MINORDER = 15;

% arguments
nargs = nargin;
if nargs < 1 || isempty( freq ), freq = 100; end
if nargs < 2 || isempty( Fs ), Fs = 1000; end
if nargs < 3 || isempty( n ), n = 0; end
if nargs < 4, band = ''; end
if nargs < 5 || isempty( tol ), tol = 0.15; end
if nargs < 6, x = []; end
if nargs < 7 || isempty( ffunc ), ffunc = 'firfilt'; end

% normalize frequencies
freq = freq / Fs * 2;
if length( freq ) ~= 2 || sum( freq < 0 ) || sum( freq == 0 ) == 2 || sum( freq > 1 )
    error( 'check FREQ/FS input' )
end

% determine filter type automatically from the freq
if isempty( band )
    if freq( 1 ) == 0
        band = 'low';
    elseif freq( 2 ) == 0
        band = 'high'; 
    else
    	band = 'bandpass'; % bandstop has to be explcitly defined
    end
end
band = lower( band );

% determine filter order automatically from the freq
nmin = ceil( FACTOR * 2 / min( freq( freq ~= 0 ) ) );
n = max( [ n nmin MINORDER ] );
odd = 1;

% design the filter
switch band
    case { 'low', 'lowpass', 'lo', 'lopass' }
        f = [ 0 freq( 2 )   freq( 2 ) * ( 1 + tol ) 1 ];
        a = [ 1 1           0                       0 ];
    case { 'high', 'highpass', 'hi', 'hipass' }
        f = [ 0 freq( 1 ) * ( 1 - tol )     freq( 1 )   1 ];
        a = [ 0 0                           1           1 ];
        odd  = 0;
    case { 'bandpass', 'band' }
        f = [ 0 freq( 1 ) * ( 1 - tol ) freq( 1 )   freq( 2 )   freq( 2 ) * ( 1 + tol ) 1   ];
        a = [ 0 0                       1           1           0                       0   ];
    case { 'stop', 'bandstop' }
        f = [ 0 freq( 1 ) * ( 1 - tol ) freq( 1 )   freq( 2 )   freq( 2 ) * ( 1 + tol ) 1   ];
        a = [ 1 1                       0           0           1                       1   ];
        odd = 0;
    otherwise
        error( 'unsuported type' )
end
if odd && ~mod( n, 2 ) || ~odd && mod( n, 2 )
    n = n + 1;
end
h = firls( n, f, a );

% actually filter
if ~isempty( x )
    if ~isa( x, 'double' )
        x = double( x );
    end
    switch lower( ffunc )
        case 'firfilt'
            xf = firfilt( x, h );
        case 'filtfilt'
            xf = zeros( size( x ) );
            for i = 1 : size( x, 2 )
                xf( :, i ) = filtfilt( h, 1, x( :, i ) );
            end
        otherwise
            xf = [];
    end
end

return

