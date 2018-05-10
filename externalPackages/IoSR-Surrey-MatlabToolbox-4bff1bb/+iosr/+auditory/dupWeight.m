function [w_itd,w_ild] = dupWeight(f)
%DUPWEIGHT Calculate duplex weighting coefficients for ITD and ILD
% 
%   [W_ITD,W_ILD] = IOSR.AUDITORY.DUPWEIGHT(F) returns the weighting
%   coefficients for ITD W_ITD and ILD W_ILD for frequency F (Hz).
% 
%   The function is derived from Ben Supper's thesis "An onset-guided
%   spatial analyser for binaural audio", with some necessary
%   modifications. Specifically, his work accounts for the ITD dominance
%   condition. This function does not, as this cannot be derived from the
%   work, since the work assumes an upper frequency limit of 13750 Hz
%   (which cannot be assumed in this function). The quadratic ramp in ITD
%   weighting is approximated linearly here.
% 
%   The input f may be an array of any size. The outputs will be the same
%   size as f, with coefficients calculated for each element.
% 
%   See also IOSR.AUDITORY.LOUDWEIGHT.

%   Copyright 2016 University of Surrey.

    %% Check input

    assert(all(f(:)>=0), 'iosr:functionalBoxPlot:invalidF', 'f should be greater than or equal to zero!')

    if any(f(:)>20000)
        warning('iosr:dupWeight:frequencyRange','Humans cannot generally hear above 20 kHz. Weighting coefficients will be set to zero.')
    end

    %% Calculate original weights using Ben's numbers

    % Frequency scale in Ben's thesis
    freq_scale = [60 150 250 350 455 570 700 845 1000 1175 ...
        1375 1600 1860 2160 2510 2925 3425 4050 4850 5850 ...
        7050 8600 10750 13750];

    % Corresponding bin numbers
    b = 1:24;

    % frequency ranges (indices)
    low = b<=8; % below cross-over region
    mid = b>8 & b<14; % cross-over region
    high = b>=14; % above cross-over region

    % pre-allocate outputs
    w_itd_orig = zeros(size(freq_scale));
    w_ild_orig = zeros(size(freq_scale));

    % Do maths
    w_itd_orig(low) = 1.597-(0.047*b(low));
    w_itd_orig(mid) = (0.0102.*(b(mid).^2))-(0.437.*b(mid))+4.06;
    w_itd_orig(high) = 0.11;
    %
    w_ild_orig(low) = 0.2;
    w_ild_orig(mid) = (0.152.*b(mid))-1.016;
    w_ild_orig(high) = 0.96;

    %% Calculate weights for input frequencies

    % ...via interpolation
    w_itd = interp1(freq_scale,w_itd_orig,f,'pchip','extrap');
    w_ild = interp1(freq_scale,w_ild_orig,f,'pchip','extrap');

    w_itd(f>20000) = 0;
    w_ild(f>20000) = 0;

end
