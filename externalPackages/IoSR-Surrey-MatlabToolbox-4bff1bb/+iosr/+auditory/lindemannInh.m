function [L_l,R_l] = lindemannInh(L,R,fs,c_inh,dim)
%LINDEMANNINH Signal pre-processing for Lindemann's cross-correlation
% 
%   [L_L,R_L] = IOSR.AUDITORY.LINDEMANNINH(L,R,FS) pre-processes left L and
%   right R signals for the cross-correlation function based on Lindemann's
%   precedence model [1,2]. A crucial parameter for Lindemann's model is
%   the "operating point", which controls the amount of inhibition. The
%   parameter is actually a function of the input related to its RMS level.
%   After half-wave rectifying and filtering the input signals L and R
%   (sampled at FS Hz) along the first non-singleton dimension, this
%   function applies a gain related to the inhibition parameter C_INH
%   (default is 0.3). The gain is identical for all rows, columns, etc.
%   Lastly, values outside of the interval [0,1] are not permitted in
%   Lindemann's model and hence these values are clipped.
% 
%   [L_L,R_L] = IOSR.AUDITORY.LINDEMANNINH(L,R,FS,C_INH) uses the specified
%   inhibition parameter C_INH. The value must be in the interval [0,1].
% 
%   [L_L,R_L] = IOSR.AUDITORY.LINDEMANNINH(L,R,FS,C_INH,DIM) pre-processes
%   L and R along the dimension DIM.
% 
%   References
%   
%   [1] Lindemann, W. (1986), Extension of a binaural cross-correlation
%       model by contralateral inhibition. I. Simulation of lateralization
%       for stationary signals, The Journal of the Acoustical Society of
%       America 80, 6, 1608-1622.
% 
%   [2] Lindemann, W. (1986), Extension of a binaural cross-correlation
%       model by contralateral inhibition. II. The law of the first wave
%       front, The Journal of the Acoustical Society of America 80, 6,
%       1623-1630.
% 
%   Further reading
% 
%   Hummersone, C., Mason, R., Brookes, T. (2013), A comparison of
%       computational precedence models for source separation in
%       reverberant environments, The Journal of the Audio Engineering
%       Society 61, 7/8, 508-520.
%   
%   See also IOSR.AUDITORY.XCORRLINDEMANN.

%   Copyright 2016 University of Surrey.

    %% check input
    
    assert(isequal(size(L),size(R)), 'iosr:lindemannInh:invalidInputs', 'L and R arrays must be the same size')
    assert(isscalar(fs), 'iosr:lindemannInh:invalidFs', 'FS must be a scalar')
    
    % default dim
    if nargin<5
        dim = find(size(L)>1,1,'first');
    end
    
    %% process
    
    % half-wave rectify
    L_l = hwr(L);
    R_l = hwr(R);

    % filter
    cutofffreq=800;
    [b,a] = butter(1,cutofffreq*2/fs);
    L_l = filter(b,a,L_l,[],dim);
    R_l = filter(b,a,R_l,[],dim);

    % inhibition parameter
    if nargin < 4
        c_inh = .3;
    else
        assert(isscalar(c_inh) & isnumeric(c_inh), 'iosr:lindemannInh:invalidCinh', 'c_inh must be a scalar');
        assert(c_inh>=0 || c_inh<=1, 'iosr:lindemannInh:invalidX', 'c_inh must be in the interval (0,1].')
    end

    % gain
    c_gamma_L = calc_gamma(L,dim);
    c_gamma_R = calc_gamma(R,dim);
    c_gamma = max([c_gamma_L(:); c_gamma_R(:)]);

    % apply parameters
    L_l = apply_gamma(L_l,c_inh,c_gamma);
    R_l = apply_gamma(R_l,c_inh,c_gamma);

    % restrict range
    L_l = clip(L_l);
    R_l = clip(R_l);

end

function y = hwr(x)
%HWR half-wave rectify

    y = max(x,0);
    
end

function gamma = calc_gamma(x,dim)
%CALC_GAMMA calculate the gamma

    gamma = sqrt(2).*iosr.dsp.rms(x,dim);
    
end

function y = apply_gamma(x,c_inh,c_gamma)
%APPLY_GAMMA apply gamma to achieve inhibition parameter

    y = (c_inh/c_gamma).*x;
    
end

function y = clip(x)
%CLIP clip input data to [0,1] interval

    y = x;
    y(y<0) = 0;
    y(y>1) = 1;

end
