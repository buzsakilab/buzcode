function rms = fastrms(x,window,dim,amp)
%FASTRMS Instantaneous root-mean-square (RMS) power via convolution.
%
%FASTRMS(X), when X is a vector, is the time-varying RMS power of X,
%computed using a 5-point rectangular window centered at each point in the
%signal. The output is the same size as X and contains, for each point in
%X, an estimate of the instantaneous power expressed in the signal.
%
%FASTRMS(X), when X is a matrix, is the time-varying RMS power of the
%columns of X. 
%
%FASTRMS(X,WINDOW), if WINDOW is a vector, computes the moving
%quadratic mean using the weights specified in WINDOW. If WINDOW is
%an integer, a LENGTH(WINDOW)-point rectangular window is used. When FASTRMS
%is being used to estimate the instantaneous amplitude of an oscillatory,
%zero-mean signal X (see below), WINDOW should be chosen based on the
%frequency content of X. Lower frequency signals require longer windows,
%whereas higher frequency signals allow shorter windows. As a rule of
%thumb, the window should be at least as long as one period of the signal. 
% 
%FASTRMS(X,WINDOW,DIM), when X is a matrix, computes the RMS power
%along the dimension DIM. (DIM specifies the "time" axis for a matrix of
%many trials.)
% 
%FASTRMS(X,WINDOW,DIM,AMP), if AMP is nonzero, applies a correction
%so that the output RMS reflects the equivalent amplitude of a sinusoidal
%input signal. That is, FASTRMS mutliplies the output by SQRT(2) to account
%for the fact that the integral of sin^2(t) over one period, t ~ [0,2*pi],
%equals (1/SQRT(2)).
% 
%The speed of FASTRMS is achieved by using convolution to compute the
%moving average of the squared signal. For this reason, FASTRMS also
%achieves maximal resolution, as the output is exactly the same size as X.
%However, the tradeoff is that some "edge effects" are incurred on the
%first and last approximately LENGTH(WINDOW)/2 samples. That is, since the
%convolution is computed using a zeropadded version of X, the RMS power
%will appear diminished near the beginning and end of the signal.
%Therefore, FASTRMS is best used on large input signals X.
% 
%EXAMPLE 
%
% Fs = 200; T = 5; N = T*Fs; t = linspace(0,T,N);
% noise = randn(N,1);
% [a,b] = butter(5, [9 12]/(Fs/2));
% x = filtfilt(a,b,noise);
% window = gausswin(0.25*Fs);
% rms = fastrms(x,window,[],1);
% plot(t,x,t,rms*[1 -1],'LineWidth',2);
% xlabel('Time (sec)'); ylabel('Signal')
% title('Instantaneous amplitude via RMS') 
% 
% 
%Created by Scott McKinney, January 2011
%http://www.mathworks.com/matlabcentral/fileexchange/authors/110216


if nargin<4
    amp=0;
end
if nargin<3 || isempty(dim)
   if size(x,1)==1 %if row vector
       dim = 2;
   else
       dim = 1; %default to computing along the columns, unless input is a row vector
   end
elseif ~isscalar(dim) || ~ismember(dim,[1 2])
    error('DIM must be the scalar value 1 or 2.')
end


if nargin < 2 || isempty(window)
    % if no window argument specified, default to length 5
    window = 5;
elseif ~isvector(window)
    error('WINDOW must be a scalar or a vector.')
end

if isscalar(window)        
    if round(window)~=window
        warning('The scalar WINDOW parameter was rounded to the nearest integer.')
    end
    window = rectwin(window);
end

if length(window) > size(x,dim);
    warning('Averaging window is longer than the length of the signal. Check inputs.')
end

if nargin < 1
    error('Input signal X is required.')
end

if length(size(x))>2
    error('FASTRMS only handles vectors and matrices.')
end

power = x.^2;

if length(size(x)) < 2 %for vectors
    rms = conv(power,window,'same');
else %for matrices
    rms = zeros(size(x));
    if dim==1
        for c = 1:size(x,2)
            rms(:,c) = conv(power(:,c),window,'same');
        end
    else
        for r = 1:size(x,1)
            rms(r,:) = conv(power(r,:),window,'same');
        end
    end                
end

rms = sqrt(rms/sum(window));
if amp
    rms = sqrt(2)*rms;
end
