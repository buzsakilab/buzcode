function [ filt_data, amp, phase ] = FiltNPhase(data,fbounds,sampfreq,varargin)
%[filt_data,amp,phase] = FiltNPhase(data,fbounds,sampfreq)
%
%Optional Arguement: filter order (num_cyc)
%
%TO DO:
%   -Add option for high/low pass filter instead of bandpass
%   -Option to input 3d data, if so, collapse into 2d plane, then re-shape output
%   -Make work with 3d matrix
%   -why problem with decimal frequencies...fix this
%   -Improve input parsing for optional inputs
%
%DLevenstein and RSwanson 2015/16
%%

%Filter order to pass to fir1. higher order improves frequency resolution
%but also increases computing time substantially.
%default: 4 cycles of low bound freq
%note: default updated to 4 from 3 cycles on 10/5/16 DL

if nargin == 4
    num_cyc = varargin{1};
else
    num_cyc = 4;
end

% Freq range possible
Hz_lowest_poss_fbound = 1/(length(data)/3/sampfreq/num_cyc);
Hz_highest_poss_found = sampfreq/2;

if fbounds(1) < Hz_lowest_poss_fbound
    error(['Lower bound of pass-band too low given number of requested cycles. Must be greater than: ' num2str(Hz_lowest_poss_fbound)])
end;
if fbounds(2) > Hz_highest_poss_found
    error(['Upper bound of pass-band exceeds Nyquist. Must be lower than: ' num2str(Hz_highest_poss_found)])
end;

%Bounds
low_bound = 1/fbounds(1);   %Low bound period - doesn't work when use decimals!!
filt_order = num_cyc*low_bound*sampfreq;    
N = ceil(filt_order);       %order must be integer
% N = 300; 

%Make bandpass FIR filter
Wn = fbounds./(0.5*sampfreq);   %Normalize frequency bounds to 1/2 sample 
                               %rate (nyquist frequency) for fir1.
B = fir1(N,Wn);                 %Designs an N'th order lowpass FIR digital
                                %filter and returns the filter coefficients
                                %in length N+1 vector B.

                                %Look at filter: fvtool(B,A)
                                
%Adjust data                                
dataROW = false;
if isrow(data)
    dataROW = true;
    data = data';
end


%Initialize
ndim = size(data,2);
filt_data = zeros(size(data)); 
amp = zeros(size(data));    
phase = zeros(size(data));    

%Filter and Hilbert transform
for dd = 1:ndim
    %Filter data and get phase using Hilbert transform
    filt_data(:,dd) = FiltFiltM(B,1,data(:,dd));
    amp(:,dd) = abs(hilbert(filt_data(:,dd)));
    phase(:,dd) = angle(hilbert(filt_data(:,dd)));
end

if dataROW
    filt_data = filt_data';
    amp = amp';
    phase = phase';
end

end
