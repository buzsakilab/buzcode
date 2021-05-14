% GAUSSKERNEL       create a 2D gaussian.
%
%                   K = GAUSSKERNEL(SIGMAX,SIGMAY,N,M)
%
%                   K                   the kernel
%                   SIGMAX, SIGMAY      std in each dimension
%                   N,M                 kernel dimensions

% 27-feb-03 ES

function K = gausskernel(sigmaX,sigmaY,N,M)
x = -(N-1)/2:(N-1)/2;
y = -(M-1)/2:(M-1)/2;
if sigmaY == 0 & sigmaX == 0
    K = zeros(M,N);
    return;
elseif sigmaY == 0
    X = -inf*ones(M,N);
    Y = X;
    X(ceil(M/2),:) = x;
    Y(ceil(M/2),:) = 0;
    sigmaY = 1;
elseif sigmaX == 0
    Y = -inf*ones(M,N);
    X = Y;
    X(:,ceil(N/2)) = 0;
    Y(:,ceil(N/2)) = y(:);
    sigmaX = 1;
else
    X = repmat(x,M,1);
    Y = repmat(y,N,1)';
end
K = 1/(2*pi*sigmaX*sigmaY)*exp(-(X.^2/2/sigmaX^2)-(Y.^2/2/sigmaY^2));
return

% to get the critical frequency of a gaussian FIR with SD of 4, use:
SD = 4; 
w = gausskernel( SD, 1, SD * 6 + 1, 1 ); w = w / sum( w );
Fs = 500; nfft = 256; pxx = abs( fft( w', nfft ) ) .^ 2; 
select = [ 1 : ( nfft + 1 ) / 2 ]; f = ( select - 1 )' * Fs / nfft; pxx = pxx( select );
s = 10*log10( pxx );
Fc = mean( f( sum( s > -3 ) + [ 0 1 ] ) )

figure, 
subplot( 1, 2, 1 ), plot( w ), 
subplot( 1, 2, 2 ), plot( f, s, '.-' ), separators( -3, [], [], 'y' )
