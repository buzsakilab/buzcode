% makegausslpfir        Gaussian LP filter
%
% win = makegausslpfir( Fc, Fs, s )
%
% Fc    corner frequency [Hz]
% Fs    sampling frequency [Hz]
% s     support [SD], minimum 3, default 4
%
% win   fir
%
% see also  firfilter

% 02-oct-13 ES

function gwin = makegausslpfir( Fc, Fs, s )

nargs = nargin;
if nargs < 1 || isempty( Fc ), Fc = 100; end
if nargs < 2 || isempty( Fs ), Fs = 1000; end
if nargs < 3 || isempty( s ), s = 4; end
s = max( s, 3 );

sd = Fs / ( 2 * pi * Fc ); 
x = -ceil( s * sd ) : ceil( s * sd ); 
gwin = 1/( 2 * pi * sd ) * exp( -( x.^2/2/sd.^2 ) ); 
gwin = gwin / sum( gwin );

return

% EOF

% to determine the proper support:
cwin = cumsum( gwin ); 
length( cwin ), [ find( cwin >= 0.001, 1, 'first' ) - 1, find( cwin <= 0.999, 1, 'last' ) + 1 ]
% 3 SD are usually enough